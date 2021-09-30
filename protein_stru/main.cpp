#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

#include "Atom.cpp"
#include "fileParser.cpp"
using std::cout;
using std::endl;
using std::string;
using std::vector;

int main(int argc, char* argv[]) {
  std::string filename = argv[1];
  cout << filename << endl;
  try {
    PDBParser parser(filename);
    parser.parse();

    parser.output2Fasta();
    parser.getL1Depth('A');
  } catch (const char* error) {
    cout << error << endl;
  }
  return 0;
}