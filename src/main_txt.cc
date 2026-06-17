//
// Build:
//   g++ -O2 -std=c++17 main_txt.cc `root-config --cflags --libs` -o repack_root_txt
//
// Run:
//   ./repack_root_txt input_list.txt output_dir
//
// Policy:
//   - input columns: Ntrig evtnum ch hitnum tdc
//   - edge is always 1
//   - output ROOT has head_TDC1 + tree_TDC1 compatible with main.cc

#include <TFile.h>
#include <TTree.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

struct TDC1Rec {
  Int_t tdc;
  Int_t edge;
  Int_t hitnum;
  Int_t evtnum;
  Int_t ch;
};

struct Head1Rec {
  Int_t gmode;
  Int_t emode;
  Int_t rmode;
  Int_t range;
  Int_t delay;
  Int_t mask[64];
  Int_t NtrigMax;
};

static void FillDefaultHead(Head1Rec& head) {
  head.gmode = 1;
  head.emode = 1;
  head.rmode = 1;
  head.range = 2;
  head.delay = 1;
  head.NtrigMax = 0;
  const int mask_values[64] = {
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0
  };
  for (int i = 0; i < 64; ++i) {
    head.mask[i] = mask_values[i];
  }
}

static int ConvertTxtToRoot(const std::string& inFile, const std::string& outDir) {
  std::filesystem::path outPath = std::filesystem::path(outDir)
                                  / std::filesystem::path(inFile).filename();
  outPath.replace_extension(".root");

  std::ifstream fin(inFile);
  if (!fin) {
    std::cerr << "[ERROR] Cannot open input file: " << inFile << "\n";
    return 2;
  }

  const std::string outFile = outPath.string();
  TFile fout(outFile.c_str(), "RECREATE");
  if (fout.IsZombie()) {
    std::cerr << "[ERROR] Cannot create output file: " << outFile << "\n";
    return 5;
  }

  Head1Rec head{};
  FillDefaultHead(head);
  TTree headTree("head_TDC1", "Head of Run - TDC1");
  headTree.Branch(
      "head1", &head,
      "gmode/I:emode/I:rmode/I:range/I:delay/I:mask[64]/I:NtrigMax/I");

  TTree dataTree("tree_TDC1", "TDC1 data");
  TDC1Rec rec{};
  dataTree.Branch("TDC1", &rec, "tdc/I:edge/I:hitnum/I:evtnum/I:ch/I");

  std::string line;
  long long lineNo = 0;
  int lastNtrig = 0;
  bool hasNtrig = false;
  while (std::getline(fin, line)) {
    ++lineNo;
    if (line.empty()) {
      continue;
    }
    std::istringstream iss(line);
    int ntrig = 0;
    int hitnum_input = 0;
    if (!(iss >> ntrig >> rec.evtnum >> rec.ch >> hitnum_input >> rec.tdc)) {
      std::cerr << "[ERROR] Failed to parse line " << lineNo << " in " << inFile
                << "\n";
      return 8;
    }
    rec.edge = 1;
    rec.hitnum = hitnum_input + 1;
    lastNtrig = ntrig;
    hasNtrig = true;
    dataTree.Fill();
  }

  if (hasNtrig) {
    head.NtrigMax = lastNtrig;
  }
  headTree.Fill();
  headTree.Write("head_TDC1");

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
    const int code = ConvertTxtToRoot(inFile, outDir);
    if (code != 0) {
      ret = code;
      break;
    }
  }

  return ret;
}
