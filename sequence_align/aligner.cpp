#include <algorithm>
#include <cassert>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <map>
using namespace std;
//https://www.youtube.com/watch?v=NqYY0PJbD3s
//https://www.youtube.com/watch?v=xyUd80lgvJc
const double inf = numeric_limits<double>::infinity();
class Cell {
 public:
  double value;
  char trace;
};

class Aligner {
 private:
  string filename1, filename2, rubricFile;
  string seqX, seqY;
  string seqXAligned, seqYAligned;
  vector<vector<int> >rubric;
  vector<vector<Cell> > M, X, Y;
  // X means seqX align to a gap in SeqY. So If choose X, it is the same as moving along SeqX (see how next in SeqX fit current SeqY)
  // Y means seqY align to a gap. So If choose Y, it is the same as moving along SeqY
  map<char, int> res2Index;
  double alignScore;
  int alignedLength, IdenticalLength;
  double sequenceIdentity;
  // M: match, X: insert @ seq1, Y: insert @ seq2
  double gapOpen, gapExt;

 public:
  Aligner(){};
  Aligner(string _filename1, string _filename2, string _rubricFile, double _gapOpen, double _gapExt)
      : filename1(_filename1),
        filename2(_filename2),
        rubricFile(_rubricFile),
        alignedLength(0),
        sequenceIdentity(0),
        IdenticalLength(0),
        gapOpen(_gapOpen),
        gapExt(_gapExt) {
    ifstream file1(_filename1), file2(_filename2);
    file1 >> seqX;
    file2 >> seqY;
    readRubric();
  };
  void readRubric();
  void align();
  double score(char a, char b);
  void print();
  void printMatrix(vector<vector<Cell> >&m);
  double calculateScore(string seqXAligned, string seqYAligned);
};
void resizeMatrix(vector<vector<Cell> >& m, int height, int width) {
  m.resize(height);
  for (int i = 0; i < height; i++) {
    m[i].resize(width);
    for (int j = 0; j < width; j++) {
      m[i][j].value = 0;
      m[i][j].trace = 'N';
      if(i==0) m[i][j].trace='X';
      if(j==0) m[i][j].trace='Y';
    }
  }
}
void Aligner::readRubric(){
  rubric = vector<vector<int> >(25, vector<int>(25));
  ifstream rubricFileReader(rubricFile);
  for(int i = 0; i < 25; i++){
    char res;
    rubricFileReader >> res;
    res2Index[res] = i;
  }
  for(int i = 0; i < 25; i++){
    for(int j = 0; j < 25; j++){
      rubricFileReader >> rubric[i][j];
    }
  }
}
void calCell(Cell& c, double m, double x, double y) {
  c.value = m;
  c.trace = 'M';
  if (x > c.value) {
    c.value = x;
    c.trace = 'X';
  }
  if (y > c.value) {
    c.value = y;
    c.trace = 'Y';
  }
}
double Aligner::score(char a, char b) {
  int indexA, indexB;
  if(res2Index.find(a)!=res2Index.end()) 
  indexA = res2Index[a];
  else indexA = res2Index['*'];
  if(res2Index.find(b)!=res2Index.end()) 
  indexB = res2Index[b];
  else indexB = res2Index['*'];
  return rubric[indexA][indexB];
}
void Aligner::align() {
  int xLength = seqX.length();
  int yLength = seqY.length();
  resizeMatrix(M, yLength + 1, xLength + 1);  // shape: yLength, xLength
  resizeMatrix(X, yLength + 1, xLength + 1);
  resizeMatrix(Y, yLength + 1, xLength + 1);
  // initialize m matrix
  M[0][0].value = 0;
  for (int i = 1; i <= yLength; i++) M[i][0].value = -inf;
  for (int j = 1; j <= xLength; j++) M[0][j].value = -inf;
  // initialize x matrix
  X[0][0].value = -gapOpen;
  // Initially for X, SeqY axis is -inf because we can't align a gap in X to a gap in Y (implied by X[0][j])
  // Moving along SeqX axis is OK because we are aligining the first j element in SeqX to gaps in SeqY
  for (int i = 1; i <= yLength; i++) X[i][0].value = -inf;
  for (int j = 1; j <= xLength; j++) X[0][j].value = X[0][j-1].value - gapExt;
  // initialize y matrix
  Y[0][0].value = -gapOpen;
  for (int i = 1; i <= yLength; i++) Y[i][0].value = Y[i-1][0].value - gapExt;
  for (int j = 1; j <= xLength; j++) Y[0][j].value =  -inf;
  // propogate weight
  for (int i = 1; i <= yLength; i++) {
    for (int j = 1; j <= xLength; j++) {
      // calculate M[i][j]
      // i-1 because the seq starts with 0
      calCell(M[i][j], M[i - 1][j - 1].value + score(seqX[j - 1], seqY[i - 1]),
              X[i - 1][j - 1].value + score(seqX[j - 1], seqY[i - 1]),
              Y[i - 1][j - 1].value + score(seqX[j - 1], seqY[i - 1]));
      // calculate X[i][j]
      calCell(X[i][j], M[i][j-1].value - gapOpen, X[i][j-1].value - gapExt,
              -inf);
      // calculate Y[i][j];
      calCell(Y[i][j], M[i-1][j].value - gapOpen, -inf,
              Y[i-1][j].value - gapExt);
    }
  }
  printMatrix(M);
  printMatrix(X);
  printMatrix(Y);
  int i = yLength, j = xLength;
  Cell curCell = M[i][j];
  char cellType = 'M';
  alignScore = M[i][j].value;
  if (X[i][j].value > curCell.value) {
    curCell = X[i][j];
    cellType = 'X';
    alignScore = X[i][j].value;
  }
  if (Y[i][j].value > curCell.value) {
    curCell = Y[i][j];
    cellType = 'Y';
    alignScore = Y[i][j].value;
  }
  cout << "start back tracing" << endl;
  while (i > 0 || j > 0) {
    if(cellType == 'X') curCell = X[i][j];
    if(cellType == 'Y') curCell = Y[i][j];
    if(cellType == 'M') curCell = M[i][j];
    if (cellType == 'M') {
      // match
      alignedLength ++;
      if(seqX[j-1]==seqY[i-1]) IdenticalLength++;
      seqXAligned += seqX[j - 1];
      seqYAligned += seqY[i - 1];
      cellType = curCell.trace;
      i--;
      j--;
      continue;
    }
    if (cellType == 'X') {
      // X aligned to a gap
      seqXAligned += seqX[j - 1];
      seqYAligned += "-";
      cellType = curCell.trace;
      j--;
      continue;
    }
    if (cellType == 'Y') {
      // Y aligned to a gap
      seqXAligned += "-";
      seqYAligned += seqY[i - 1];
      cellType = curCell.trace;
      i--;
      continue;
    }
  }
  reverse(seqXAligned.begin(), seqXAligned.end());
  reverse(seqYAligned.begin(), seqYAligned.end());
}
void Aligner::print() {
  cout << "Aligned Sequence Score:" << alignScore << endl;
  cout << "Aligned Length: "<< alignedLength << endl;
  cout << "Identical Length:" << IdenticalLength << endl;
  cout << "Sequence identity: " << static_cast<double>(IdenticalLength)/static_cast<double>(seqY.length()) << endl; 
  cout << "Calculated Score:" << calculateScore(seqXAligned, seqYAligned) << endl;
  cout << "Aligned Sequence Result:" << endl;
  int resPerLine = 50;
  assert(seqXAligned.length() == seqYAligned.length());
  int start = 0;
  while (start + resPerLine < seqXAligned.length()) {
    cout << "> X: " << seqXAligned.substr(start, start + resPerLine) << endl;
    cout << "> Y: " << seqYAligned.substr(start, start + resPerLine) << endl;
    cout << endl;
    start += resPerLine;
  }
  cout << "> X: " << seqXAligned.substr(start, seqXAligned.length() - start)
       << endl;
  cout << "> Y: " << seqYAligned.substr(start, seqYAligned.length() - start)
       << endl;
}
void Aligner::printMatrix(vector<vector<Cell> >&m){
    cout << "Start printMatrix: \n";
    for(int i = 0; i < m.size(); i++){
        for(int j = 0; j < m[0].size(); j++){
            cout << std::setw(5) << m[i][j].value << ' ';
        }
        cout << endl;
    }
    cout << "End printMatrix: \n";

}
double Aligner::calculateScore(string seqXAligned, string seqYAligned){
  int seqLength = seqXAligned.length();
  bool isGapOpenX = true, isGapOpenY = true;
  double alignScore = 0; 
  for(int i = 0; i < seqLength; i++){
    // cout << "aligned Score =" << alignScore << endl;
    if(seqXAligned[i] == '-' || seqYAligned[i] == '-'){
      if(seqXAligned[i] == '-' && isGapOpenX){
        alignScore -= gapOpen;
        isGapOpenX = false;
      }
      else if(seqYAligned[i] == '-' && isGapOpenY){
        alignScore -= gapOpen;
        isGapOpenY = false;
      }
      else alignScore -= gapExt;
    }
    else{
      alignScore += score(seqXAligned[i], seqYAligned[i]);
      isGapOpenX = true;
      isGapOpenY = true;
    }
  }
  return alignScore;
}
