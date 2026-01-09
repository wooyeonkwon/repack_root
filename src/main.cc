//
// Build:
//   g++ -O2 -std=c++17 main.cc repack_root.cc `root-config --cflags --libs` -o repack_root
//
// Run:
//   ./repack_root input.root output.root
//

#include "repack_root.h"

int main(int argc, char** argv) {
  return RunRepackRoot(argc, argv);
}
