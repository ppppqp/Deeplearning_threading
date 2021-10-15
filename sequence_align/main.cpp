#include <iostream>

#include "aligner.cpp"
int main() {
  Aligner aligner("seqX", "seqY", 3, 1);
  aligner.align();
  aligner.print();
  return 0;
}
using namespace std;