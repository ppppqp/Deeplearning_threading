#include <iostream>

#include "aligner.cpp"
int main() {
  Aligner aligner("seqX", "seqY", "BLOSUM62",  11, 1);
  cout << "calculating...\n";
  aligner.align();
  aligner.print();
  return 0;
}
using namespace std;