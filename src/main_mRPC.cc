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
//   - first leading hits locate each channel peak with 2000-ps bins
//   - a local +/-3000-ps histogram with 200-ps bins is used for the Gaussian fit
//   - invalid fits warn and use a local median/MAD window, or the fixed peak-search window
//   - fallback windows stay inside the local +/-3000-ps peak-search region
//   - TimingCalibration records the method and actual selection bounds per channel
//   - each raw hit stores its channel timing center in offset and calibrated raw time in tdc_cali_raw
//   - clean hits additionally store their timing center in offset_clean and tdc = tdc_raw - offset
//   - evtnum and *_raw branches are always written, even when nHit==0
//   - streaming grouping by evtnum (assumes evtnum monotonic in file)

#include <TFile.h>
#include <TF1.h>
#include <TFitResult.h>
#include <TFitResultPtr.h>
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
  // raw input values, with per-raw-hit calibration values aligned by index
  std::vector<int> ch_raw, tdc_raw, edge_raw, hitnum_raw;
  std::vector<double> offset, tdc_cali_raw;

  // clean (hitnum==1 && edge==1), with mRPC channel mapping applied
  std::vector<int> ch;
  std::vector<double> tdc, offset_clean;
  std::vector<int> clean_rawIndex;

  // number of hits after clean selection and threshold removal
  int nHit = 0;

  void reset() {
    ch_raw.clear();
    tdc_raw.clear();
    edge_raw.clear();
    hitnum_raw.clear();
    offset.clear();
    tdc_cali_raw.clear();
    ch.clear();
    tdc.clear();
    offset_clean.clear();
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
  std::vector<double> windowScalePs;
  std::vector<double> cleanLower;
  std::vector<double> cleanUpper;
  // 0: no first leading hits; 1: Gaussian; 2: local median/MAD; 3: fixed window.
  std::vector<int> method;
};

static double Median(std::vector<double> values) {
  if (values.empty()) return std::numeric_limits<double>::quiet_NaN();
  std::sort(values.begin(), values.end());
  const size_t mid = values.size() / 2;
  return values.size() % 2 ? values[mid] : 0.5 * (values[mid - 1] + values[mid]);
}

static ChannelCalibration ComputeChannelCalibration(TTree* tin, TDC1Rec& rec) {
  const int kMaxCh = 64;
  const double kBinWidth = 2000.0;
  const double kFitBinWidth = 200.0;
  const double kInitialSigma = 1000.0;
  const double kFitHalfWidth = 3.0 * kInitialSigma;

  ChannelCalibration calibration;
  calibration.offsets.assign(kMaxCh + 1, 0.0);
  calibration.windowScalePs.assign(kMaxCh + 1, 0.0);
  calibration.method.assign(kMaxCh + 1, 0);
  calibration.cleanLower.assign(kMaxCh + 1, std::numeric_limits<double>::quiet_NaN());
  calibration.cleanUpper.assign(kMaxCh + 1, std::numeric_limits<double>::quiet_NaN());

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
    if (rec.hitnum != 1 || rec.edge != 1) {
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
    hTdcRaw.SetDirectory(nullptr);
    for (const int tdcRaw : tdcs) {
      hTdcRaw.Fill(tdcRaw);
    }

    const double modeTdc = hTdcRaw.GetBinCenter(hTdcRaw.GetMaximumBin());
    const double fitMin = modeTdc - kFitHalfWidth;
    const double fitMax = modeTdc + kFitHalfWidth;
    const int nFitBins = static_cast<int>(std::lround((fitMax - fitMin) / kFitBinWidth));
    const std::string fineName = histName + "_peak";
    TH1D hPeak(fineName.c_str(), "Channel timing peak;TDC time (ps);Counts",
               nFitBins, fitMin, fitMax);
    hPeak.SetDirectory(nullptr);
    int peakEntries = 0;
    std::vector<double> peakTimes;
    for (const int tdcRaw : tdcs) {
      if (tdcRaw >= fitMin && tdcRaw < fitMax) {
        hPeak.Fill(tdcRaw);
        ++peakEntries;
        peakTimes.push_back(tdcRaw);
      }
    }

    TF1 gausFit(fitName.c_str(), "gaus", fitMin, fitMax);
    gausFit.SetParameters(hPeak.GetMaximum(), modeTdc, kInitialSigma);
    gausFit.SetParLimits(1, fitMin, fitMax);
    gausFit.SetParLimits(2, 1.0, kFitHalfWidth);
    // Poisson likelihood includes empty bins; I averages the model within each bin.
    // S retains diagnostics, N avoids storing/drawing the temporary function.
    TFitResultPtr result = hPeak.Fit(&gausFit, "QSNRLI");
    const int status = result;
    const double mean = gausFit.GetParameter(1);
    const double sigma = gausFit.GetParameter(2);
    const double boundaryTolerance = 0.1; // ps
    const bool valid = result.Get() && status == 0 && result->IsValid()
        && result->CovMatrixStatus() == 3 && result->Ndf() > 0
        && std::isfinite(gausFit.GetParameter(0)) && gausFit.GetParameter(0) > 0.0
        && std::isfinite(mean) && mean > fitMin + boundaryTolerance
        && mean < fitMax - boundaryTolerance
        && std::isfinite(sigma) && sigma > 1.0 + boundaryTolerance
        && sigma < kFitHalfWidth - boundaryTolerance;
    std::cout << "[CALIB] ch=" << ch << " firstLeadingHits=" << tdcs.size()
              << " peakEntries=" << peakEntries << " bins=" << nFitBins
              << " status=" << status
              << " covStatus=" << (result.Get() ? result->CovMatrixStatus() : -1)
              << " mean_ps=" << mean << " sigma_ps=" << sigma
              << " accepted=" << valid << "\n";
    if (!valid) {
      // Estimate only the local peak, never the long raw-time tail. This is a
      // selection fallback, not a Gaussian resolution measurement. Require at
      // least 10 local hits for the robust estimate; otherwise keep the original
      // fixed search window. Both choices remain explicitly flagged in output.
      const double median = Median(peakTimes);
      std::vector<double> deviations;
      for (const double time : peakTimes) deviations.push_back(std::abs(time - median));
      const double robustScale = 1.4826 * Median(deviations);
      const bool useRobust = peakTimes.size() >= 10 && std::isfinite(median)
          && std::isfinite(robustScale) && robustScale >= kFitBinWidth / 2.0
          && robustScale < kFitHalfWidth;
      const double center = useRobust ? median : modeTdc;
      const double scale = useRobust ? robustScale : kInitialSigma;
      calibration.method[ch] = useRobust ? 2 : 3;
      calibration.offsets[ch] = center;
      calibration.windowScalePs[ch] = scale;
      calibration.cleanLower[ch] = std::max(fitMin, center - 3.0 * scale);
      calibration.cleanUpper[ch] = std::min(fitMax, center + 3.0 * scale);
      std::cerr << "[WARN] Invalid timing fit for channel " << ch
                << "; using " << (useRobust ? "local median/MAD" : "fixed peak window")
                << ": offset_ps=" << center << " window_scale_ps=" << scale
                << " lower_ps=" << calibration.cleanLower[ch]
                << " upper_ps=" << calibration.cleanUpper[ch]
                << ". Fallback scale is NOT a measured Gaussian sigma.\n";
      continue;
    }

    calibration.offsets[ch] = mean;
    calibration.windowScalePs[ch] = sigma;
    calibration.method[ch] = 1;
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

  const ChannelCalibration calibration = ComputeChannelCalibration(tin, rec);
  std::cout << "[INFO] mRPC per-channel timing calibration computed\n";

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
  // Persist fallback information so it survives beyond the console log.
  TTree calibrationTree("TimingCalibration", "Methods: 0=empty, 1=Gaussian, 2=local median/MAD, 3=fixed window");
  int calCh = 0, calMethod = 0;
  double calOffset = 0.0, calScale = 0.0, calSigma = 0.0, calLower = 0.0, calUpper = 0.0;
  calibrationTree.Branch("ch", &calCh);
  calibrationTree.Branch("method", &calMethod);
  calibrationTree.Branch("offset_ps", &calOffset);
  calibrationTree.Branch("window_scale_ps", &calScale);
  calibrationTree.Branch("gaussian_sigma_ps", &calSigma);
  calibrationTree.Branch("lower_ps", &calLower);
  calibrationTree.Branch("upper_ps", &calUpper);
  for (calCh = 1; calCh <= 64; ++calCh) {
    calMethod = calibration.method[calCh];
    calOffset = calibration.offsets[calCh];
    calScale = calibration.windowScalePs[calCh];
    calSigma = calMethod == 1 ? calScale : std::numeric_limits<double>::quiet_NaN();
    calLower = calibration.cleanLower[calCh];
    calUpper = calibration.cleanUpper[calCh];
    calibrationTree.Fill();
  }
  TTree tout("Events", "Event-level repacked TDC1 for mRPC");

  int o_evtnum = 0;
  EventBuffers ev;

  tout.Branch("evtnum", &o_evtnum);

  tout.Branch("ch_raw", &ev.ch_raw);
  tout.Branch("tdc_raw", &ev.tdc_raw);
  tout.Branch("edge_raw", &ev.edge_raw);
  tout.Branch("hitnum_raw", &ev.hitnum_raw);
  tout.Branch("offset", &ev.offset);
  tout.Branch("tdc_cali_raw", &ev.tdc_cali_raw);

  tout.Branch("ch", &ev.ch);
  tout.Branch("tdc", &ev.tdc);
  tout.Branch("offset_clean", &ev.offset_clean);
  tout.Branch("clean_rawIndex", &ev.clean_rawIndex);

  tout.Branch("nHit", &ev.nHit);

  ev.reset();
  int curEvt = -1;
  Long64_t firstLeadingHits = 0;
  Long64_t acceptedHits = 0;

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

    const int mappedCh = MapMrpcChannel(rec.ch);
    const bool knownCh = mappedCh >= 1 && mappedCh < static_cast<int>(calibration.offsets.size());
    const double offset = knownCh ? calibration.offsets[mappedCh] : 0.0;
    ev.offset.push_back(offset);
    ev.tdc_cali_raw.push_back(static_cast<double>(rec.tdc) - offset);

    if (knownCh && calibration.method[mappedCh] != 0
        && rec.hitnum == 1 && rec.edge == 1) {
      ++firstLeadingHits;
      const double cleanLower = calibration.cleanLower[mappedCh];
      const double cleanUpper = calibration.cleanUpper[mappedCh];
      if (static_cast<double>(rec.tdc) >= cleanLower
          && static_cast<double>(rec.tdc) <= cleanUpper) {
        ++acceptedHits;
        ev.ch.push_back(mappedCh);
        ev.tdc.push_back(static_cast<double>(rec.tdc) - offset);
        ev.offset_clean.push_back(offset);
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

  std::cout << "[SELECTION] rawRecords=" << nEnt
            << " firstLeadingHits=" << firstLeadingHits
            << " acceptedHits=" << acceptedHits << "\n";
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
