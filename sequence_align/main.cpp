#include <iostream>

#include "aligner.cpp"
int main() {
  Aligner aligner("seqX", "seqY", "BLOSUM62",  12, 1);
  cout << "calculating...\n";
  aligner.align();
  aligner.print();
  string refX = "EQA-HFALFFNQGQCCCAGSRTFVQEDIYAEFVERSVARAKSRVVGN-----PFDSRTEQGPQVDETQFKKVLGYI";
  string refY = "KSGKEEGLKLLCGGGAAADRGYFIQPTVFGD-LQDGMTIAKEEIFGPVMQILKFKSMEEVVGRANNSKYG-----L";
  cout << "ref score is: : "<< aligner.calculateScore(refX, refY) << endl;
  return 0;
}
using namespace std;