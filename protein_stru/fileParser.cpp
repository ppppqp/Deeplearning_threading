#ifndef parser_cpp
#define parser_cpp
#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "Atom.cpp"
#include "Chain.cpp"
#include "Residue.cpp"
#include "utils.hpp"
using std::ceil;
using std::cout;
using std::map;
using std::string;
using std::vector;
const double RESIDUALS_PER_LINE = 13;
string dropSemiColon(string s) {
  if (s[s.length() - 1] == ';') s.erase(s.length() - 1, 1);
  return s;
}
class firstLineSwitch {
 public:
  bool title;
  bool caveat;
  bool compnd;
  bool source;
  firstLineSwitch() : title(true), caveat(true), compnd(true), source(true) {}
};
class PDBParser {
 private:
  char uniformChainNum;
  string filename;
  string proteinName;
  map<char, Chain> chains;  // chains it contains
  vector<char> chainIdList;
  firstLineSwitch Switch;
  std::ifstream fin;
  // molecule -> chain -> residue -> atom
 public:
  PDBParser(string _filename)
      : filename(_filename){};
  void parse();
  void parseCompnd(std::stringstream& ss);
  void parseAtom(string line);
  void parseSeqres(std::stringstream& ss);
  void parseHeader(std::stringstream& ss);
  void output2Fasta();
  void output2PDB(char chainNum, string prefix);
  void inferChain(char chainNum, int residueId);
  void inferResidue(char chainNum, int residueId, string type);
  void getL1Depth(char chainNum);
};

void PDBParser::parse() {
  fin.open(filename.c_str());
  if (!fin.is_open()) {
    throw "Fail to open file";
  }
  string line;
  string junk;
  int count = 0;
  size_t nameStart = filename.find_last_of('/');
  proteinName = filename.substr(nameStart + 1, 4);
  std::transform(proteinName.begin(), proteinName.end(), proteinName.begin(),
                ::tolower);
  if(filename.size() > 8) uniformChainNum = filename[nameStart + 5];
  while (getline(fin, line)) {  // successful load file
    std::stringstream ss(line);
    string lineType;
    ss >> lineType;
    // if (lineType == "COMPND") {
    //   parseCompnd(ss);
    //   continue;
    // }
    if (lineType == "SEQRES") {  // parse seqres field
      // parseSeqres(ss);
      continue;
    }
    count++;
    if (lineType == "HEADER") {
      // parseHeader(ss);
      continue;
    }
    if (lineType == "ATOM") {  // parse ATOM field
      parseAtom(line);
      continue;
    }
  }
  map<char, Chain>::iterator it;
  for(it = chains.begin(); it!=chains.end(); it++) it->second.cleanNonCARes();
};
void PDBParser::parseHeader(std::stringstream& ss) {
  while (ss >> proteinName) {
  };
  std::transform(proteinName.begin(), proteinName.end(), proteinName.begin(),
                 ::tolower);
}
void PDBParser::parseCompnd(std::stringstream& ss) {
  int molId;
  string line, lineType, junk;
  if (Switch.compnd) {  // parse molId
    ss >> lineType >> junk >> molId;
    Switch.compnd = false;
  } else {
    ss >> lineType >> junk >> junk >> molId;
  }
  getline(fin, line);  // parse molecule name
  string moleculeName = line.substr(line.find(':') + 2, line.find(';'));
  getline(fin, line);  // parse chain name
  string chainName = line.substr(line.find(':') + 2, line.find(';'));
  getline(fin, line);  // parse synonym name
  string synonym = line.substr(line.find(':') + 2, line.find(';'));
  getline(fin, line);  // parse engineered
  string engineered = line.substr(line.find(':') + 2, line.find(';'));

  vector<string> chainNum;
  std::stringstream is(chainName);
  string chain;
  while (getline(is, chain, ',')) chainNum.push_back(chain);
  for (int i = 0; i < chainNum.size(); i++) {
    Chain c;
    c.id = molId;
    c.name = dropSemiColon(moleculeName);
    c.synonym = dropSemiColon(synonym);
    c.engineered = dropSemiColon(engineered);
    if (chainNum[i].size() == 1)
      c.chainNum = chainNum[i][0];
    else
      c.chainNum = chainNum[i][1];
    chains[c.chainNum] = c;
    chainIdList.push_back(c.chainNum);
  }
}
void PDBParser::parseAtom(string line) {
  AtomInfo atomInfo;
  bool debug = false;
  if(debug) cout << line << endl;
  atomInfo.atomType = line.substr(12,4);
  trim(atomInfo.atomType);

  atomInfo.residue = line.substr(17,3);
  trim(atomInfo.residue);

  atomInfo.chainNum = line[21];
  if(atomInfo.chainNum== ' ') atomInfo.chainNum = uniformChainNum;
  
  atomInfo.residueId =  stoi(line.substr(22,4));

  atomInfo.x = stod(line.substr(30,8));
  atomInfo.y = stod(line.substr(38, 8));
  atomInfo.z = stod(line.substr(46,8));
  string eleStr = "";
  if(line.length() > 60) atomInfo.occupancy = stod(line.substr(54,6));
  if(line.length() > 66) atomInfo.beta = stod(line.substr(60,6));
  if(line.length() >= 78 && line.substr(76,2)!="  ") {
      eleStr = line.substr(76,2);
      trim(eleStr);
      if(eleStr.length()!=0) atomInfo.element = eleStr[0];
  }
  else {
      if(!(atomInfo.atomType[0]>='0' && atomInfo.atomType[0] <= '9'))
          atomInfo.element = atomInfo.atomType[0];
      else atomInfo.element = atomInfo.atomType[1];
  }
  Atom atom(atomInfo);
  if(debug) cout << "break point 1" << endl;
  inferChain(atomInfo.chainNum, atomInfo.residueId);  // make sure chain exists
  atomInfo.residueId -= chains[atomInfo.chainNum].resBaseIndex; // convert to zero based index
  inferResidue(atomInfo.chainNum, atomInfo.residueId,
               atomInfo.residue);  // make sure residue exists
  if(debug) cout << "break point 3" << endl;
  chains[atomInfo.chainNum].residues[atomInfo.residueId].atoms.push_back(
      atom);
}
void PDBParser::inferChain(char chainNum, int residueId) {
  map<char, Chain>::iterator it;
  it = chains.find(chainNum);
  if (it == chains.end()) {
    Chain c(chainNum, residueId);
    chainIdList.push_back(chainNum);
    chains[chainNum] = c;
  }
  return;
}
void PDBParser::inferResidue(char chainNum, int residueId, string type) {
  vector<Residue>& resVec = chains[chainNum].residues;
  if (resVec.size() < residueId + 1) {  // residueId starts from 0
    // for example, if current size is 6, residueId is 7, then we have X, X, X, X,
    // X, X, _, X so we need to push one dummy item. After we push the dummy
    // item, size=7, residueId=7
    int gap = residueId - resVec.size();
    for (int i = 0; i < gap; i++) {
      // push until residueId == resVec.size()+1
      Residue r(false);
      resVec.push_back(r);
    }
    // residueId == resVec.size()
    Residue r(type, residueId, chainNum);
    resVec.push_back(r);
    chains[chainNum].validResNum++;
  } else if (!resVec[residueId].valid) {  // in case disorder
    Residue r(type, residueId, chainNum);
    resVec[residueId] = r;
    chains[chainNum].validResNum++;
  }
}
void PDBParser::output2PDB(char chainNum, string prefix){
  std::ofstream fout;

  string fileName = prefix + proteinName + chainNum;
  cout << "filename is " << fileName << endl;
  fout.open(fileName.c_str());
  Chain & chain = chains[chainNum];
  vector<Residue> &residues = chain.residues;
  int lineNum = 1;
  for(int i = 0; i < residues.size(); i++){
    if(residues[i].valid){
      vector<Atom> &atoms = residues[i].atoms;
      for(int j = 0; j < atoms.size(); j++){
        // cout << "i = " << i << " j = " << j << endl;
        if(atoms[j].valid){
          Atom& atom = atoms[j];
          fout << "ATOM  "; 
          fout << setw(5) << lineNum << " ";
          if(atom.atomType.size() < 4) fout << " " << std::left << setw(3) << atom.atomType;
          else fout << atom.atomType;
          fout << std::right << setw(4) << atom.residue;
          
          fout << setw(4) << atom.residueId << "   ";
          fout << setw(8) << std::fixed << std::setprecision(3) << atom.x;
          fout << setw(8) << std::fixed << std::setprecision(3) << atom.y;
          fout << setw(8) << std::fixed << std::setprecision(3) << atom.z;

          fout << setw(6) << std::fixed << std::setprecision(2) <<atom.occupancy;
          fout << setw(6) << std::fixed << std::setprecision(2) << atom.beta;
          fout << "          ";
          fout << "  " << atom.element << endl;
          // fout << "  " << atomInfo.charge << endl;
          lineNum++;
        }
      }
    }
  }
}
void PDBParser::output2Fasta() {
  std::ofstream fout;

  const int resPerLine = 70;
  fout.open("pdb2fasta");
  for (int i = 0; i < chainIdList.size(); i++) {
    int count = 1;
    Chain& c = chains[chainIdList[i]];
    fout << '>' << proteinName << chainIdList[i] << "  " << c.validResNum
         << '\n';
    for (int j = 0; j < c.residues.size(); j++) {
      if (c.residues[j].valid) {
        fout << pdb2fasta(c.residues[j].type);
        if (count % resPerLine == 0) fout << '\n';
        count++;
      }
    }
    fout << '\n';
  }
  fout.close();
}
void PDBParser::parseSeqres(std::stringstream& ss) {
  int lineNum;
  char chainNum;
  string line, lineType;
  int ResidualNum;
  ss >> lineNum >> chainNum >> ResidualNum;
  int totalLineNum =
      ceil(static_cast<double>(ResidualNum) / RESIDUALS_PER_LINE);
  vector<string> lineSeq;
  lineSeq.reserve(ResidualNum);
  string res;
  int count = 1;  // one count per chain
  while (ss >> res) {
    lineSeq.push_back(res);  // first line
    Residue r;
    r.residueId = count;
    r.chainId = chainNum;
    r.type = res;
    chains[chainNum].residues.push_back(r);
    count++;
  }
  for (int i = 1; i < totalLineNum; i++) {
    getline(fin, line);  // rest line
    std::stringstream ss2(line);
    ss2 >> lineType >> lineNum >> chainNum >> ResidualNum;
    while (ss2 >> res) {
      lineSeq.push_back(res);
      Residue r;
      r.residueId = count;
      r.chainId = chainNum;
      r.type = res;
      chains[chainNum].residues.push_back(r);
      count++;
    }
  }
  chains[chainNum].print();
}
void PDBParser::getL1Depth(char chainNum) {
  chains[chainNum].L1Depth("residue-residue", "global");
}
#endif