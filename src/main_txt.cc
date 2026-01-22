//
// Build:
//   g++ -O2 -std=c++17 main_txt.cc `root-config --cflags --libs` -o repack_root_txt
//
// Run:
//   ./repack_root_txt input_list.txt output_dir
//
// Policy:
//   - clean = (hitnum==1 && edge==1)  [leading only]
//   - isPaired computed on clean hits, opp = 65 - ch
//   - left: 1..32, right: 33..64
//   - streaming grouping by evtnum (assumes evtnum monotonic in file)

#include <TFile.h>
#include <TTree.h>

#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <limits>
#include <iostream>
#include <string>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <sstream>

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
  std::vector<int> clean_rawIndex;
  std::vector<unsigned char> isPaired;

  // flags
  bool hasLeft = false;
  bool hasRight = false;
  unsigned char trigCategory = 0;   // 0 none, 1 left, 2 right, 3 both
  bool isCoincidenceEvent = false;

  // fastest
  int fast_ch = -1;
  int fast_tdc = -1;
  int fast_cleanIndex = -1;
  int fast_rawIndex = -1;
  bool fast_isPaired = false;

  // QA
  bool hasMultiHit = false;
  int nMultiHit = 0;

  void reset() {
    ch_raw.clear(); tdc_raw.clear(); edge_raw.clear(); hitnum_raw.clear();
    ch.clear(); tdc.clear(); clean_rawIndex.clear(); isPaired.clear();

    hasLeft = false; hasRight = false;
    trigCategory = 0;
    isCoincidenceEvent = false;

    fast_ch = -1; fast_tdc = -1;
    fast_cleanIndex = -1; fast_rawIndex = -1;
    fast_isPaired = false;

    hasMultiHit = false; nMultiHit = 0;
  }
};

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
  if (!pairedCh.empty()) {
    for (size_t i = 0; i < ev.ch.size(); ++i) {
      if (pairedCh.find(ev.ch[i]) != pairedCh.end()) ev.isPaired[i] = 1;
    }
  }
  ev.isCoincidenceEvent = !pairedCh.empty();

  // fastest (min tdc among CLEAN)
  int bestT = std::numeric_limits<int>::max();
  int bestIdx = -1;
  for (size_t i = 0; i < ev.tdc.size(); ++i) {
    if (ev.tdc[i] < bestT) {
      bestT = ev.tdc[i];
      bestIdx = static_cast<int>(i);
    }
  }

  if (bestIdx >= 0) {
    ev.fast_cleanIndex = bestIdx;
    ev.fast_ch = ev.ch[bestIdx];
    ev.fast_tdc = ev.tdc[bestIdx];
    ev.fast_rawIndex = ev.clean_rawIndex[bestIdx];
    ev.fast_isPaired = (ev.isPaired[bestIdx] != 0);
  } else {
    ev.fast_cleanIndex = -1;
    ev.fast_ch = -1;
    ev.fast_tdc = -1;
    ev.fast_rawIndex = -1;
    ev.fast_isPaired = false;
  }
}

static int ProcessFile(const std::string& inFile, const std::string& outDir) {
  std::filesystem::path outPath = std::filesystem::path(outDir)
                                  / std::filesystem::path(inFile).filename();
  outPath.replace_extension(".root");

  std::ifstream fin(inFile);
  if (!fin) {
    std::cerr << "[ERROR] Cannot open input file: " << inFile << "\n";
    return 2;
  }

  // Output
  const std::string outFile = outPath.string();
  TFile fout(outFile.c_str(), "RECREATE");
  if (fout.IsZombie()) {
    std::cerr << "[ERROR] Cannot create output file: " << outFile << "\n";
    return 5;
  }

  int ntrig_max = 0;
  TTree head("head_TDC1", "head_TDC1");
  head.Branch("NtrigMax", &ntrig_max);

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
  tout.Branch("clean_rawIndex", &ev.clean_rawIndex);
  tout.Branch("isPaired", &ev.isPaired);

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

  // Streaming loop
  std::unordered_map<int, EventBuffers> events;
  events.reserve(4096);
  std::string line;
  long long lineNo = 0;
  while (std::getline(fin, line)) {
    ++lineNo;
    if (line.empty()) {
      continue;
    }
    std::istringstream iss(line);
    TDC1Rec rec{};
    int ntrig = 0;
    if (!(iss >> ntrig >> rec.evtnum >> rec.ch >> rec.hitnum >> rec.tdc)) {
      std::cerr << "[ERROR] Failed to parse line " << lineNo << " in " << inFile
                << "\n";
      return 8;
    }
    rec.edge = 1;
    if (ntrig > ntrig_max) {
      ntrig_max = ntrig;
    }
    EventBuffers& curEvent = events[rec.evtnum];

    // push raw
    const int rawIndex = static_cast<int>(curEvent.ch_raw.size());
    curEvent.ch_raw.push_back(rec.ch);
    curEvent.tdc_raw.push_back(rec.tdc);
    curEvent.edge_raw.push_back(rec.edge);
    curEvent.hitnum_raw.push_back(rec.hitnum);

    // QA
    if (rec.hitnum > 1) {
      curEvent.hasMultiHit = true;
      curEvent.nMultiHit++;
    }

    // push clean (hitnum==1 && edge==1)
    if (rec.hitnum == 1 && rec.edge == 1) {
      curEvent.ch.push_back(rec.ch);
      curEvent.tdc.push_back(rec.tdc);
      curEvent.clean_rawIndex.push_back(rawIndex);
    }
  }

  std::vector<int> evtnums;
  evtnums.reserve(events.size());
  for (const auto& entry : events) {
    evtnums.push_back(entry.first);
  }
  std::sort(evtnums.begin(), evtnums.end());

  for (int evtnum : evtnums) {
    EventBuffers& curEvent = events[evtnum];
    FinalizeEvent(evtnum, curEvent);
    o_evtnum = evtnum;
    tout.Fill();
  }

  head.Fill();
  head.Write("head_TDC1");
  fout.Write();
  fout.Close();

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
