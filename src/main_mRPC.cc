//
// Build:
//   g++ -O2 -std=c++17 main_mRPC.cc `root-config --cflags --libs` -o repack_root_mRPC
//
// Run:
//   ./repack_root_mRPC input_list.txt output_dir
//
// Policy:
//   - clean = (hitnum==1 && edge==1)  [leading only]
//   - mRPC reads only one side, so no pair/coincidence branches are produced
//   - input channel 1..32 is reversed for output clean channel: ch = 33 - ch_raw
//   - tdc_raw is binned in 1000-unit bins and fit with a Gaussian for each channel
//   - each channel keeps clean hits inside that channel's Gaussian mean +/- 3 sigma window
//   - each channel's Gaussian mean is stored as offset and tdc = tdc_raw - offset
//   - evtnum and *_raw branches are always written, even when nHit==0
//   - streaming grouping by evtnum (assumes evtnum monotonic in file)

#include <TFile.h>
#include <TF1.h>
#include <TH1D.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

struct TDC1Rec {
  Int_t tdc;
  Int_t edge;
  Int_t hitnum;
  Int_t evtnum;
  Int_t ch;
};

struct EventBuffers {
  // raw input values
  std::vector<int> ch_raw, tdc_raw, edge_raw, hitnum_raw;

  // clean (hitnum==1 && edge==1), with mRPC channel mapping applied
  std::vector<int> ch;
  std::vector<double> tdc, offset;
  std::vector<int> clean_rawIndex;

  // number of hits after clean selection and threshold removal
  int nHit = 0;

  void reset() {
    ch_raw.clear();
    tdc_raw.clear();
    edge_raw.clear();
    hitnum_raw.clear();
    ch.clear();
    tdc.clear();
    offset.clear();
    clean_rawIndex.clear();
    nHit = 0;
  }
};

static long long CountEvtnumDrops(TTree* tin, TDC1Rec& rec) {
  const Long64_t n = tin->GetEntries();
  int prev = -1;
  long long drops = 0;
  for (Long64_t i = 0; i < n; ++i) {
    tin->GetEntry(i);
    if (prev != -1 && rec.evtnum < prev) drops++;
    prev = rec.evtnum;
  }
  return drops;
}

static inline int MapMrpcChannel(int chRaw) {
  if (chRaw >= 1 && chRaw <= 32) {
    return 33 - chRaw;
  }
  return chRaw;
}

struct ChannelCalibration {
  std::vector<double> offsets;
  std::vector<double> gaussianSigma;
  std::vector<double> cleanLower;
  std::vector<double> cleanUpper;
};

static ChannelCalibration ComputeChannelCalibration(TTree* tin, TDC1Rec& rec) {
  const int kMaxCh = 64;
  const double kBinWidth = 20000.0;

  ChannelCalibration calibration;
  calibration.offsets.assign(kMaxCh + 1, 0.0);
  calibration.gaussianSigma.assign(kMaxCh + 1, 0.0);
  calibration.cleanLower.assign(kMaxCh + 1, std::numeric_limits<double>::lowest());
  calibration.cleanUpper.assign(kMaxCh + 1, std::numeric_limits<double>::max());

  const Long64_t n = tin->GetEntries();
  if (n <= 0) {
    return calibration;
  }

  std::vector<std::vector<int>> tdcRawByCh(kMaxCh + 1);
  for (Long64_t i = 0; i < n; ++i) {
    tin->GetEntry(i);
    if (rec.ch < 1 || rec.ch > kMaxCh) {
      continue;
    }
    tdcRawByCh[MapMrpcChannel(rec.ch)].push_back(rec.tdc);
  }

  for (int ch = 1; ch <= kMaxCh; ++ch) {
    const std::vector<int>& tdcs = tdcRawByCh[ch];
    if (tdcs.empty()) {
      continue;
    }

    const auto [minIt, maxIt] = std::minmax_element(tdcs.begin(), tdcs.end());
    const double histMin = std::floor(static_cast<double>(*minIt) / kBinWidth) * kBinWidth;
    const double histMax = (std::floor(static_cast<double>(*maxIt) / kBinWidth) + 1.0) * kBinWidth;
    const int nBins = std::max(1, static_cast<int>((histMax - histMin) / kBinWidth));

    const std::string histName = "h_mrpc_tdc_raw_ch" + std::to_string(ch);
    const std::string fitName = "mrpc_tdc_raw_gaus_ch" + std::to_string(ch);
    TH1D hTdcRaw(histName.c_str(), "mRPC tdc_raw;tdc_raw;counts", nBins, histMin, histMax);
    for (const int tdcRaw : tdcs) {
      hTdcRaw.Fill(tdcRaw);
    }

    double mean = hTdcRaw.GetMean();
    double sigma = hTdcRaw.GetRMS();

    if (hTdcRaw.GetEntries() > 0) {
      const int modeBin = hTdcRaw.GetMaximumBin();
      const double modeEntryCount = hTdcRaw.GetBinContent(modeBin);
      const double modeTdc = hTdcRaw.GetBinCenter(modeBin);

      TF1 gausFit(fitName.c_str(), "gaus", histMin, histMax);
      gausFit.SetParameters(modeEntryCount, modeTdc, kBinWidth);
      hTdcRaw.Fit(&gausFit, "Q0");

      const double fitSigma = std::abs(gausFit.GetParameter(2));
      if (fitSigma > 0.0) {
        mean = gausFit.GetParameter(1);
        sigma = fitSigma;
      }
    }

    calibration.offsets[ch] = mean;
    calibration.gaussianSigma[ch] = sigma;
    if (sigma > 0.0) {
      calibration.cleanLower[ch] = mean - 3.0 * sigma;
      calibration.cleanUpper[ch] = mean + 3.0 * sigma;
    }
  }

  return calibration;
}

static int ProcessFile(const std::string& inFile, const std::string& outDir) {
  const std::filesystem::path outPath = std::filesystem::path(outDir)
                                        / std::filesystem::path(inFile).filename();

  TFile fin(inFile.c_str(), "READ");
  if (fin.IsZombie()) {
    std::cerr << "[ERROR] Cannot open input file: " << inFile << "\n";
    return 2;
  }

  TTree* tin = (TTree*)fin.Get("tree_TDC1");
  if (!tin) {
    std::cerr << "[ERROR] tree_TDC1 not found\n";
    return 3;
  }

  TDC1Rec rec;
  if (tin->SetBranchAddress("TDC1", &rec) < 0) {
    std::cerr << "[ERROR] SetBranchAddress(\"TDC1\") failed\n";
    return 4;
  }

  const long long drops = CountEvtnumDrops(tin, rec);
  if (drops > 0) {
    std::cerr << "[WARN] evtnum is not monotonic (drops=" << drops
              << "). Streaming grouping may break.\n";
  }

  const std::string outFile = outPath.string();
  TFile fout(outFile.c_str(), "RECREATE");
  if (fout.IsZombie()) {
    std::cerr << "[ERROR] Cannot create output file: " << outFile << "\n";
    return 5;
  }

  if (TTree* head = (TTree*)fin.Get("head_TDC1")) {
    fout.cd();
    TTree* headOut = head->CloneTree(-1, "fast");
    headOut->Write("head_TDC1");
  }

  fout.cd();
  TTree tout("Events", "Event-level repacked TDC1 for mRPC");

  int o_evtnum = 0;
  EventBuffers ev;

  tout.Branch("evtnum", &o_evtnum);

  tout.Branch("ch_raw", &ev.ch_raw);
  tout.Branch("tdc_raw", &ev.tdc_raw);
  tout.Branch("edge_raw", &ev.edge_raw);
  tout.Branch("hitnum_raw", &ev.hitnum_raw);

  tout.Branch("ch", &ev.ch);
  tout.Branch("tdc", &ev.tdc);
  tout.Branch("offset", &ev.offset);
  tout.Branch("clean_rawIndex", &ev.clean_rawIndex);

  tout.Branch("nHit", &ev.nHit);

  const ChannelCalibration calibration = ComputeChannelCalibration(tin, rec);
  std::cout << "[INFO] mRPC per-channel Gaussian calibration computed\n";

  ev.reset();
  int curEvt = -1;

  const Long64_t nEnt = tin->GetEntries();
  for (Long64_t i = 0; i < nEnt; ++i) {
    tin->GetEntry(i);

    if (curEvt == -1) curEvt = rec.evtnum;

    if (rec.evtnum != curEvt) {
      ev.nHit = static_cast<int>(ev.ch.size());
      o_evtnum = curEvt;
      tout.Fill();
      ev.reset();
      curEvt = rec.evtnum;
    }

    const int rawIndex = static_cast<int>(ev.ch_raw.size());
    ev.ch_raw.push_back(rec.ch);
    ev.tdc_raw.push_back(rec.tdc);
    ev.edge_raw.push_back(rec.edge);
    ev.hitnum_raw.push_back(rec.hitnum);

    if (rec.hitnum == 1 && rec.edge == 1) {
      const int mappedCh = MapMrpcChannel(rec.ch);
      const bool knownCh = mappedCh >= 1 && mappedCh < static_cast<int>(calibration.offsets.size());
      const double offset = knownCh ? calibration.offsets[mappedCh] : 0.0;
      const double cleanLower = knownCh ? calibration.cleanLower[mappedCh]
                                       : std::numeric_limits<double>::lowest();
      const double cleanUpper = knownCh ? calibration.cleanUpper[mappedCh]
                                       : std::numeric_limits<double>::max();
      if (static_cast<double>(rec.tdc) >= cleanLower
          && static_cast<double>(rec.tdc) <= cleanUpper) {
        ev.ch.push_back(mappedCh);
        ev.tdc.push_back(static_cast<double>(rec.tdc) - offset);
        ev.offset.push_back(offset);
        ev.clean_rawIndex.push_back(rawIndex);
      }
    }
  }

  if (curEvt != -1) {
    ev.nHit = static_cast<int>(ev.ch.size());
    o_evtnum = curEvt;
    tout.Fill();
  }

  fout.Write();
  fout.Close();
  fin.Close();

  std::cout << "[INFO] Done. Output: " << outFile << "\n";
  return 0;
}

int main(int argc, char** argv) {
  if (argc < 3) {
    std::cerr << "Usage: " << argv[0] << " input_list.txt output_dir\n";
    return 1;
  }

  const std::string listFile = argv[1];
  const std::string outDir = argv[2];

  std::filesystem::path outDirPath(outDir);
  if (!std::filesystem::exists(outDirPath)) {
    std::error_code ec;
    if (!std::filesystem::create_directories(outDirPath, ec)) {
      std::cerr << "[ERROR] Cannot create output directory: " << outDir
                << " (" << ec.message() << ")\n";
      return 6;
    }
  }

  std::ifstream inList(listFile);
  if (!inList) {
    std::cerr << "[ERROR] Cannot open input list file: " << listFile << "\n";
    return 7;
  }

  std::string inFile;
  int ret = 0;
  while (std::getline(inList, inFile)) {
    if (inFile.empty()) {
      continue;
    }
    const int code = ProcessFile(inFile, outDir);
    if (code != 0) {
      ret = code;
      break;
    }
  }

  return ret;
}
