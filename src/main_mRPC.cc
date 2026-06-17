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
//   - channel offsets and an upper clean threshold are estimated from single-hit TDC values
//   - clean hits beyond the threshold are removed; events with nHit==0 are not written
//   - streaming grouping by evtnum (assumes evtnum monotonic in file)

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <unordered_map>
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
  double cleanThreshold = std::numeric_limits<double>::max();
};

static ChannelCalibration ComputeChannelCalibration(TTree* tin, TDC1Rec& rec) {
  const int kMaxCh = 32;
  const int kModeBinWidth = 1000;
  const int kTargetHitCount = 1;

  struct CleanHit {
    int ch;
    int tdc;
  };

  std::vector<std::vector<int>> eventTdcs;
  std::vector<int> allTdcs;
  eventTdcs.reserve(1024);
  allTdcs.reserve(4096);

  std::vector<CleanHit> eventHits;
  eventHits.reserve(kMaxCh);
  int curEvt = -1;

  auto flushThresholdEvent = [&]() {
    std::vector<int> tdcs;
    tdcs.reserve(eventHits.size());
    for (const CleanHit& h : eventHits) {
      tdcs.push_back(h.tdc);
      allTdcs.push_back(h.tdc);
    }
    eventTdcs.push_back(std::move(tdcs));
  };

  const Long64_t n = tin->GetEntries();
  for (Long64_t i = 0; i < n; ++i) {
    tin->GetEntry(i);
    if (curEvt == -1) {
      curEvt = rec.evtnum;
    }
    if (rec.evtnum != curEvt) {
      flushThresholdEvent();
      eventHits.clear();
      curEvt = rec.evtnum;
    }
    if (rec.hitnum == 1 && rec.edge == 1 && rec.ch >= 1 && rec.ch <= kMaxCh) {
      eventHits.push_back({MapMrpcChannel(rec.ch), rec.tdc});
    }
  }
  if (curEvt != -1) {
    flushThresholdEvent();
  }

  int selectedThreshold = std::numeric_limits<int>::max();
  if (!allTdcs.empty()) {
    std::sort(allTdcs.begin(), allTdcs.end());
    allTdcs.erase(std::unique(allTdcs.begin(), allTdcs.end()), allTdcs.end());
    std::sort(allTdcs.rbegin(), allTdcs.rend());

    long long bestEventCount = -1;
    int bestThreshold = allTdcs.front();
    for (const int threshold : allTdcs) {
      long long eventCountWithTarget = 0;
      for (const auto& tdcs : eventTdcs) {
        int matched = 0;
        for (const int v : tdcs) {
          if (v < threshold) {
            matched++;
          }
        }
        if (matched == kTargetHitCount) {
          eventCountWithTarget++;
        }
      }
      if (eventCountWithTarget > bestEventCount
          || (eventCountWithTarget == bestEventCount && threshold > bestThreshold)) {
        bestEventCount = eventCountWithTarget;
        bestThreshold = threshold;
      }
    }
    selectedThreshold = bestThreshold;
  }

  std::vector<int> minVal(kMaxCh + 1, std::numeric_limits<int>::max());
  std::vector<std::unordered_map<int, int>> binCounts(kMaxCh + 1);

  auto isTargetEvent = [&](const std::vector<CleanHit>& hits) {
    int matched = 0;
    for (const CleanHit& h : hits) {
      if (h.tdc < selectedThreshold) {
        matched++;
      }
    }
    return matched == kTargetHitCount;
  };

  auto accumulateEventForMode = [&](const std::vector<CleanHit>& hits) {
    if (!isTargetEvent(hits)) {
      return;
    }
    for (const CleanHit& h : hits) {
      if (h.tdc >= selectedThreshold) {
        continue;
      }
      minVal[h.ch] = std::min(minVal[h.ch], h.tdc);
      const int binIdx = h.tdc / kModeBinWidth;
      binCounts[h.ch][binIdx]++;
    }
  };

  eventHits.clear();
  curEvt = -1;
  for (Long64_t i = 0; i < n; ++i) {
    tin->GetEntry(i);
    if (curEvt == -1) {
      curEvt = rec.evtnum;
    }
    if (rec.evtnum != curEvt) {
      accumulateEventForMode(eventHits);
      eventHits.clear();
      curEvt = rec.evtnum;
    }
    if (rec.hitnum == 1 && rec.edge == 1 && rec.ch >= 1 && rec.ch <= kMaxCh) {
      eventHits.push_back({MapMrpcChannel(rec.ch), rec.tdc});
    }
  }
  if (curEvt != -1) {
    accumulateEventForMode(eventHits);
  }

  std::vector<double> modeVal(kMaxCh + 1, 0.0);
  for (int ch = 1; ch <= kMaxCh; ++ch) {
    int bestCount = 0;
    int bestBinIdx = 0;
    for (const auto& kv : binCounts[ch]) {
      if (kv.second > bestCount || (kv.second == bestCount && kv.first < bestBinIdx)) {
        bestCount = kv.second;
        bestBinIdx = kv.first;
      }
    }
    if (bestCount > 0) {
      modeVal[ch] = static_cast<double>(bestBinIdx * kModeBinWidth + (kModeBinWidth / 2));
    }
  }

  std::vector<long long> sumVal(kMaxCh + 1, 0);
  std::vector<long long> countVal(kMaxCh + 1, 0);
  auto accumulateEventForAverage = [&](const std::vector<CleanHit>& hits) {
    if (!isTargetEvent(hits)) {
      return;
    }
    for (const CleanHit& h : hits) {
      if (h.tdc >= selectedThreshold || binCounts[h.ch].empty()) {
        continue;
      }
      const int minTdc = minVal[h.ch];
      const double modeTdc = modeVal[h.ch];
      const double maxTdc = 2.0 * modeTdc - static_cast<double>(minTdc);
      if (static_cast<double>(h.tdc) >= static_cast<double>(minTdc)
          && static_cast<double>(h.tdc) <= maxTdc) {
        sumVal[h.ch] += h.tdc;
        countVal[h.ch] += 1;
      }
    }
  };

  eventHits.clear();
  curEvt = -1;
  for (Long64_t i = 0; i < n; ++i) {
    tin->GetEntry(i);
    if (curEvt == -1) {
      curEvt = rec.evtnum;
    }
    if (rec.evtnum != curEvt) {
      accumulateEventForAverage(eventHits);
      eventHits.clear();
      curEvt = rec.evtnum;
    }
    if (rec.hitnum == 1 && rec.edge == 1 && rec.ch >= 1 && rec.ch <= kMaxCh) {
      eventHits.push_back({MapMrpcChannel(rec.ch), rec.tdc});
    }
  }
  if (curEvt != -1) {
    accumulateEventForAverage(eventHits);
  }

  ChannelCalibration calibration;
  calibration.offsets.assign(kMaxCh + 1, 0.0);
  calibration.cleanThreshold = static_cast<double>(selectedThreshold);
  for (int ch = 1; ch <= kMaxCh; ++ch) {
    if (countVal[ch] > 0) {
      calibration.offsets[ch] = static_cast<double>(sumVal[ch]) / static_cast<double>(countVal[ch]);
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
  std::cout << "[INFO] mRPC clean threshold: " << calibration.cleanThreshold << "\n";

  ev.reset();
  int curEvt = -1;

  const Long64_t nEnt = tin->GetEntries();
  for (Long64_t i = 0; i < nEnt; ++i) {
    tin->GetEntry(i);

    if (curEvt == -1) curEvt = rec.evtnum;

    if (rec.evtnum != curEvt) {
      ev.nHit = static_cast<int>(ev.ch.size());
      if (ev.nHit > 0) {
        o_evtnum = curEvt;
        tout.Fill();
      }
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
      const double offset = (mappedCh >= 1 && mappedCh < static_cast<int>(calibration.offsets.size()))
        ? calibration.offsets[mappedCh]
        : 0.0;
      if (static_cast<double>(rec.tdc) < calibration.cleanThreshold) {
        ev.ch.push_back(mappedCh);
        ev.tdc.push_back(rec.tdc);
        ev.offset.push_back(offset);
        ev.clean_rawIndex.push_back(rawIndex);
      }
    }
  }

  if (curEvt != -1) {
    ev.nHit = static_cast<int>(ev.ch.size());
    if (ev.nHit > 0) {
      o_evtnum = curEvt;
      tout.Fill();
    }
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
