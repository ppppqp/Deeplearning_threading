#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
// ./output ../../PDB/9rubA.pdb
#include "Atom.cpp"
#include "fileParser.cpp"
using std::cout;
using std::endl;
using std::string;
using std::vector;

int main(int argc, char* argv[]) {
  std::string filename = argv[1];
  bool isChain = false;
  if(argc > 1) isChain = true;
  cout << filename << endl;
  try {
    PDBParser parser(filename, isChain);
    parser.parse();

    parser.output2Fasta();
    parser.getL1Depth('A');
  } catch (const char* error) {
    cout << error << endl;
  }
  return 0;
}