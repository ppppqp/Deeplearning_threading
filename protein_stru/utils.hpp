#ifndef Utils_hpp
#define Utils_hpp
#include <iostream>
#include <string>
using std::cout;
using std::endl;
using std::string;
#include <algorithm>
#include <cctype>
#include <string>
// ssh -Y your_username@lighthouse.arc-ts.umich.edu
// HHpred PfamScan
char pdb2fasta(string pbdResidue) {
  if (pbdResidue == "ALA") return 'A';
  if (pbdResidue == "ARG") return 'R';
  if (pbdResidue == "ASN") return 'N';
  if (pbdResidue == "ASP") return 'D';
  if (pbdResidue == "CYS") return 'C';
  if (pbdResidue == "GLU") return 'E';
  if (pbdResidue == "GLN") return 'Q';
  if (pbdResidue == "GLY") return 'G';
  if (pbdResidue == "HIS") return 'H';
  if (pbdResidue == "HYP") return 'O';
  if (pbdResidue == "IIE") return 'I';
  if (pbdResidue == "LEU") return 'L';
  if (pbdResidue == "LYS") return 'K';
  if (pbdResidue == "MET") return 'M';
  if (pbdResidue == "PHE") return 'F';
  if (pbdResidue == "PRO") return 'P';
  if (pbdResidue == "SER") return 'S';
  if (pbdResidue == "THR") return 'T';
  if (pbdResidue == "TRP") return 'W';
  if (pbdResidue == "TYR") return 'Y';
  if (pbdResidue == "VAL") return 'V';
  if (pbdResidue == "ILE") return 'I';
  if (pbdResidue == "DT") return 'T';
  if (pbdResidue == "DG") return 'G';
  if (pbdResidue == "DA") return 'A';
  if (pbdResidue == "DC") return 'C';

  return 'G';
//   cout << pbdResidue << endl;
//   throw "No corresponding pdb2fasta conversion";
}

// trim from start (in place)
static inline void ltrim(std::string &s) {
  s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](unsigned char ch) {
            return !std::isspace(ch);
          }));
}

// trim from end (in place)
static inline void rtrim(std::string &s) {
  s.erase(std::find_if(s.rbegin(), s.rend(),
                       [](unsigned char ch) { return !std::isspace(ch); })
              .base(),
          s.end());
}

// trim from both ends (in place)
static inline void trim(std::string &s) {
  ltrim(s);
  rtrim(s);
}

#endif