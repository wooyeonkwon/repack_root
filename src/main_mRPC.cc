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
//   - streaming grouping by evtnum (assumes evtnum monotonic in file)

#include <TFile.h>
#include <TTree.h>

#include <filesystem>
#include <fstream>
#include <iostream>
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
  std::vector<int> ch, tdc;
  std::vector<double> offset;
  std::vector<int> clean_rawIndex;

  // QA
  bool hasMultiHit = false;
  int nMultiHit = 0;

  void reset() {
    ch_raw.clear();
    tdc_raw.clear();
    edge_raw.clear();
    hitnum_raw.clear();
    ch.clear();
    tdc.clear();
    offset.clear();
    clean_rawIndex.clear();
    hasMultiHit = false;
    nMultiHit = 0;
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

  tout.Branch("hasMultiHit", &ev.hasMultiHit);
  tout.Branch("nMultiHit", &ev.nMultiHit);

  ev.reset();
  int curEvt = -1;

  const Long64_t nEnt = tin->GetEntries();
  for (Long64_t i = 0; i < nEnt; ++i) {
    tin->GetEntry(i);

    if (curEvt == -1) curEvt = rec.evtnum;

    if (rec.evtnum != curEvt) {
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

    if (rec.hitnum > 1) {
      ev.hasMultiHit = true;
      ev.nMultiHit++;
    }

    if (rec.hitnum == 1 && rec.edge == 1) {
      ev.ch.push_back(MapMrpcChannel(rec.ch));
      ev.tdc.push_back(rec.tdc);
      ev.offset.push_back(0.0);
      ev.clean_rawIndex.push_back(rawIndex);
    }
  }

  if (curEvt != -1) {
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
