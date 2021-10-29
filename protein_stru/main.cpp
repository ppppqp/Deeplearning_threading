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
////./output ../data/pdb/2PU3.pdb
int main(int argc, char* argv[]) {
  std::string filename = argv[1];
  cout << filename << endl;
  try {
    PDBParser parser(filename);
    parser.parse();
    parser.output2Fasta();
    parser.output2PDB('A', "./");
    parser.getL1Depth('A');
  } catch (const char* error) {
    cout << error << endl;
  }
  return 0;
}