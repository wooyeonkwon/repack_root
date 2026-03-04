//
// Build:
//   g++ -O2 -std=c++17 main.cc `root-config --cflags --libs` -o repack_root
//
// Run:
//   ./repack_root input.root output.root
//
// Policy:
//   - clean = (hitnum==1 && edge==1)  [leading only]
//   - isPaired computed on clean hits, opp = 65 - ch
//   - left: 1..32, right: 33..64
//   - streaming grouping by evtnum (assumes evtnum monotonic in file)

#include <TFile.h>
#include <TTree.h>

#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <limits>
#include <iostream>
#include <string>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <cmath>

struct TDC1Rec {
  Int_t tdc;
  Int_t edge;
  Int_t hitnum;
  Int_t evtnum;
  Int_t ch;
};

static inline bool IsLeft(int ch)  { return (ch >= 1 && ch <= 32); }
static inline bool IsRight(int ch) { return (ch >= 33 && ch <= 64); }
static inline int  OppChannel(int ch) { return 65 - ch; }

struct EventBuffers {
  // raw
  std::vector<int> ch_raw, tdc_raw, edge_raw, hitnum_raw;

  // clean (hitnum==1 && edge==1)
  std::vector<int> ch, tdc;
  std::vector<double> offset, tdc_cali;
  std::vector<int> clean_rawIndex;
  std::vector<unsigned char> isPaired;
  std::vector<int> pairIdx;

  // flags
  bool hasLeft = false;
  bool hasRight = false;
  unsigned char trigCategory = 0;   // 0 none, 1 left, 2 right, 3 both
  bool isCoincidenceEvent = false;

  // fastest
  int fast_ch = -1;
  double fast_tdc = -1.0;
  int fast_cleanIndex = -1;
  int fast_rawIndex = -1;
  bool fast_isPaired = false;

  // QA
  bool hasMultiHit = false;
  int nMultiHit = 0;

  void reset() {
    ch_raw.clear(); tdc_raw.clear(); edge_raw.clear(); hitnum_raw.clear();
    ch.clear(); tdc.clear(); offset.clear(); tdc_cali.clear();
    clean_rawIndex.clear(); isPaired.clear(); pairIdx.clear();

    hasLeft = false; hasRight = false;
    trigCategory = 0;
    isCoincidenceEvent = false;

    fast_ch = -1; fast_tdc = -1.0;
    fast_cleanIndex = -1; fast_rawIndex = -1;
    fast_isPaired = false;

    hasMultiHit = false; nMultiHit = 0;
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

static std::vector<double> ComputeChannelOffsets(TTree* tin, TDC1Rec& rec) {
  const int kMaxCh = 64;
  const int kModeBinWidth = 1000;
  const int kPairHitCountTarget = 2;
  const int kPairCountTarget = kPairHitCountTarget / 2;

  struct CleanHit {
    int ch;
    int tdc;
  };

  struct PairSumData {
    std::vector<int> sums;
    std::unordered_map<int, int> chToPairSum;
  };

  auto BuildPairSumData = [&](const std::vector<CleanHit>& hits) {
    PairSumData out;
    std::unordered_map<int, int> chToTdc;
    chToTdc.reserve(hits.size());
    for (const CleanHit& h : hits) {
      chToTdc.emplace(h.ch, h.tdc);
    }

    for (const auto& kv : chToTdc) {
      const int ch = kv.first;
      const int opp = OppChannel(ch);
      if (ch >= opp) {
        continue;
      }
      auto itOpp = chToTdc.find(opp);
      if (itOpp == chToTdc.end()) {
        continue;
      }
      const int tdcSum = kv.second + itOpp->second;
      out.sums.push_back(tdcSum);
      out.chToPairSum[ch] = tdcSum;
      out.chToPairSum[opp] = tdcSum;
    }
    return out;
  };

  std::vector<std::vector<int>> eventPairSums;
  std::vector<int> allPairSums;
  eventPairSums.reserve(1024);
  allPairSums.reserve(4096);

  std::vector<CleanHit> eventHits;
  eventHits.reserve(64);
  int curEvt = -1;

  const Long64_t n = tin->GetEntries();
  for (Long64_t i = 0; i < n; ++i) {
    tin->GetEntry(i);
    if (curEvt == -1) {
      curEvt = rec.evtnum;
    }
    if (rec.evtnum != curEvt) {
      PairSumData pairData = BuildPairSumData(eventHits);
      eventPairSums.push_back(pairData.sums);
      allPairSums.insert(allPairSums.end(), pairData.sums.begin(), pairData.sums.end());
      eventHits.clear();
      curEvt = rec.evtnum;
    }
    if (rec.hitnum == 1 && rec.edge == 1 && rec.ch >= 1 && rec.ch <= kMaxCh) {
      eventHits.push_back({rec.ch, rec.tdc});
    }
  }
  if (curEvt != -1) {
    PairSumData pairData = BuildPairSumData(eventHits);
    eventPairSums.push_back(pairData.sums);
    allPairSums.insert(allPairSums.end(), pairData.sums.begin(), pairData.sums.end());
  }

  int selectedThreshold = std::numeric_limits<int>::max();
  if (!allPairSums.empty()) {
    std::sort(allPairSums.begin(), allPairSums.end());
    allPairSums.erase(std::unique(allPairSums.begin(), allPairSums.end()), allPairSums.end());
    std::sort(allPairSums.rbegin(), allPairSums.rend());

    long long bestEventCount = -1;
    int bestThreshold = allPairSums.front();
    for (const int threshold : allPairSums) {
      long long eventCountWithTarget = 0;
      for (const auto& sums : eventPairSums) {
        int matched = 0;
        for (const int v : sums) {
          if (v < threshold) {
            matched++;
          }
        }
        if (matched == kPairCountTarget) {
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

  auto IsTargetEvent = [&](const PairSumData& pairData) {
    int matched = 0;
    for (const int v : pairData.sums) {
      if (v < selectedThreshold) {
        matched++;
      }
    }
    // matched counts pair-sums (one value per pair), so target is 1 pair (=2 hits).
    return matched == kPairCountTarget;
  };

  auto accumulateEventForMode = [&](const std::vector<CleanHit>& hits) {
    if (hits.empty()) {
      return;
    }
    PairSumData pairData = BuildPairSumData(hits);
    if (!IsTargetEvent(pairData)) {
      return;
    }
    for (const CleanHit& h : hits) {
      auto it = pairData.chToPairSum.find(h.ch);
      if (it == pairData.chToPairSum.end()) {
        continue;
      }
      if (it->second >= selectedThreshold) {
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
      eventHits.push_back({rec.ch, rec.tdc});
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
    if (hits.empty()) {
      return;
    }
    PairSumData pairData = BuildPairSumData(hits);
    if (!IsTargetEvent(pairData)) {
      return;
    }
    for (const CleanHit& h : hits) {
      auto it = pairData.chToPairSum.find(h.ch);
      if (it == pairData.chToPairSum.end()) {
        continue;
      }
      if (it->second >= selectedThreshold) {
        continue;
      }
      if (binCounts[h.ch].empty()) {
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
      eventHits.push_back({rec.ch, rec.tdc});
    }
  }
  if (curEvt != -1) {
    accumulateEventForAverage(eventHits);
  }

  std::vector<double> offsets(kMaxCh + 1, 0.0);
  for (int ch = 1; ch <= kMaxCh; ++ch) {
    if (countVal[ch] > 0) {
      offsets[ch] = static_cast<double>(sumVal[ch]) / static_cast<double>(countVal[ch]);
    }
  }
  return offsets;
}

static void FinalizeEvent(int evtnum, EventBuffers& ev) {
  // hasLeft/hasRight from CLEAN
  ev.hasLeft = false;
  ev.hasRight = false;
  for (int c : ev.ch) {
    if (IsLeft(c))  ev.hasLeft = true;
    if (IsRight(c)) ev.hasRight = true;
  }
  ev.trigCategory = static_cast<unsigned char>((ev.hasLeft ? 1 : 0) + (ev.hasRight ? 2 : 0));

  // paired on CLEAN
  std::unordered_set<int> chset;
  chset.reserve(ev.ch.size());
  for (int c : ev.ch) chset.insert(c);

  std::unordered_set<int> pairedCh;
  pairedCh.reserve(chset.size());
  for (int c : chset) {
    const int opp = OppChannel(c);
    if (chset.find(opp) != chset.end()) {
      pairedCh.insert(c);
      pairedCh.insert(opp);
    }
  }

  ev.isPaired.assign(ev.ch.size(), 0);
  ev.pairIdx.assign(ev.ch.size(), -1);
  std::unordered_map<int, int> chToCleanIdx;
  chToCleanIdx.reserve(ev.ch.size());
  for (size_t i = 0; i < ev.ch.size(); ++i) {
    chToCleanIdx.emplace(ev.ch[i], static_cast<int>(i));
  }
  if (!pairedCh.empty()) {
    for (size_t i = 0; i < ev.ch.size(); ++i) {
      if (pairedCh.find(ev.ch[i]) != pairedCh.end()) {
        ev.isPaired[i] = 1;
        const int opp = OppChannel(ev.ch[i]);
        auto itPair = chToCleanIdx.find(opp);
        if (itPair != chToCleanIdx.end()) {
          ev.pairIdx[i] = itPair->second;
        }
      }
    }
  }
  ev.isCoincidenceEvent = !pairedCh.empty();

  // fastest (min calibrated tdc among CLEAN)
  double bestT = std::numeric_limits<double>::max();
  int bestIdx = -1;
  for (size_t i = 0; i < ev.tdc_cali.size(); ++i) {
    if (ev.tdc_cali[i] < bestT) {
      bestT = ev.tdc_cali[i];
      bestIdx = static_cast<int>(i);
    }
  }

  if (bestIdx >= 0) {
    ev.fast_cleanIndex = bestIdx;
    ev.fast_ch = ev.ch[bestIdx];
    ev.fast_tdc = ev.tdc_cali[bestIdx];
    ev.fast_rawIndex = ev.clean_rawIndex[bestIdx];
    ev.fast_isPaired = (ev.isPaired[bestIdx] != 0);
  } else {
    ev.fast_cleanIndex = -1;
    ev.fast_ch = -1;
    ev.fast_tdc = -1.0;
    ev.fast_rawIndex = -1;
    ev.fast_isPaired = false;
  }
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

  // Check monotonicity (informational)
  const long long drops = CountEvtnumDrops(tin, rec);
  if (drops > 0) {
    std::cerr << "[WARN] evtnum is not monotonic (drops=" << drops
              << "). Streaming grouping may break.\n";
  }

  // Output
  const std::string outFile = outPath.string();
  TFile fout(outFile.c_str(), "RECREATE");
  if (fout.IsZombie()) {
    std::cerr << "[ERROR] Cannot create output file: " << outFile << "\n";
    return 5;
  }

  // Clone header tree if exists
  if (TTree* head = (TTree*)fin.Get("head_TDC1")) {
    fout.cd();
    TTree* headOut = head->CloneTree(-1, "fast");
    headOut->Write("head_TDC1");
  }

  // Events Tree
  fout.cd();
  TTree tout("Events", "Event-level repacked TDC1 (raw + clean leading)");

  int o_evtnum = 0;
  EventBuffers ev;

  // Branches
  tout.Branch("evtnum", &o_evtnum);

  tout.Branch("ch_raw", &ev.ch_raw);
  tout.Branch("tdc_raw", &ev.tdc_raw);
  tout.Branch("edge_raw", &ev.edge_raw);
  tout.Branch("hitnum_raw", &ev.hitnum_raw);

  tout.Branch("ch", &ev.ch);
  tout.Branch("tdc", &ev.tdc);
  tout.Branch("offset", &ev.offset);
  tout.Branch("tdc_cali", &ev.tdc_cali);
  tout.Branch("clean_rawIndex", &ev.clean_rawIndex);
  tout.Branch("isPaired", &ev.isPaired);
  tout.Branch("pairIdx", &ev.pairIdx);

  tout.Branch("hasLeft", &ev.hasLeft);
  tout.Branch("hasRight", &ev.hasRight);
  tout.Branch("trigCategory", &ev.trigCategory);
  tout.Branch("isCoincidenceEvent", &ev.isCoincidenceEvent);

  tout.Branch("fast_ch", &ev.fast_ch);
  tout.Branch("fast_tdc", &ev.fast_tdc);
  tout.Branch("fast_cleanIndex", &ev.fast_cleanIndex);
  tout.Branch("fast_rawIndex", &ev.fast_rawIndex);
  tout.Branch("fast_isPaired", &ev.fast_isPaired);

  tout.Branch("hasMultiHit", &ev.hasMultiHit);
  tout.Branch("nMultiHit", &ev.nMultiHit);

  const std::vector<double> channelOffsets = ComputeChannelOffsets(tin, rec);

  // Streaming loop
  ev.reset();
  int curEvt = -1;

  const Long64_t nEnt = tin->GetEntries();
  for (Long64_t i = 0; i < nEnt; ++i) {
    tin->GetEntry(i);

    if (curEvt == -1) curEvt = rec.evtnum;

    // event boundary
    if (rec.evtnum != curEvt) {
      // finalize + fill
      FinalizeEvent(curEvt, ev);
      o_evtnum = curEvt;
      tout.Fill();

      // reset for next
      ev.reset();
      curEvt = rec.evtnum;
    }

    // push raw
    const int rawIndex = static_cast<int>(ev.ch_raw.size());
    ev.ch_raw.push_back(rec.ch);
    ev.tdc_raw.push_back(rec.tdc);
    ev.edge_raw.push_back(rec.edge);
    ev.hitnum_raw.push_back(rec.hitnum);

    // QA
    if (rec.hitnum > 1) {
      ev.hasMultiHit = true;
      ev.nMultiHit++;
    }

    // push clean (hitnum==1 && edge==1)
    if (rec.hitnum == 1 && rec.edge == 1) {
      ev.ch.push_back(rec.ch);
      ev.tdc.push_back(rec.tdc);
      const double offset = (rec.ch >= 1 && rec.ch < static_cast<int>(channelOffsets.size()))
        ? channelOffsets[rec.ch]
        : 0.0;
      ev.offset.push_back(offset);
      ev.tdc_cali.push_back(static_cast<double>(rec.tdc) - offset);
      ev.clean_rawIndex.push_back(rawIndex);
    }
  }

  // last event
  if (curEvt != -1) {
    FinalizeEvent(curEvt, ev);
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
