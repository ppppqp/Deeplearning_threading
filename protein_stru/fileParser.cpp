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
  vector<Residue>altRes;
  int lastResId;
  std::ifstream fin;
  // molecule -> chain -> residue -> atom
 public:
  PDBParser(string _filename)
      : filename(_filename), lastResId(-1){};
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
  void cleanAltRes();
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
      // cout << count << ' ' << line << endl;
      parseAtom(line);
      continue;
    }
    
  }
  map<char, Chain>::iterator it;
  cleanAltRes();
  for(it = chains.begin(); it!=chains.end(); it++) it->second.cleanNonCARes();
};
void PDBParser::cleanAltRes(){
  for(int i = 0; i < altRes.size(); i++){
      if(altRes[i].getAtom("CA").valid){
        if(altRes[i].type.length() > 3){
          altRes[i].type.erase(0,1);// get rid of the alt tag
        }
        inferResidue(altRes[i].chainId, altRes[i].residueId, altRes[i].type);  // make sure residue exists
        chains[altRes[i].chainId].residues[altRes[i].residueId].atoms = altRes[i].atoms;
        break;
      }
    }
    altRes.clear();
}
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
  char altTag = 'X';

  atomInfo.atomType = line.substr(12,4);
  trim(atomInfo.atomType);
  if(atomInfo.atomType.length() == 4 && (atomInfo.atomType[0] == 'A' || atomInfo.atomType[0] == 'B')){
    altTag = atomInfo.atomType[0];
    atomInfo.atomType.erase(0,1);
  }

  atomInfo.residue = line.substr(16,4);
  trim(atomInfo.residue);


  atomInfo.chainNum = line[21];
  if(atomInfo.chainNum== ' ') atomInfo.chainNum = uniformChainNum;
  
  string residueIdStr = line.substr(22,4);
  for(int i = 0; i < residueIdStr.length(); i++){
    if(isalpha(residueIdStr[i])){
      altTag = residueIdStr[i];
      residueIdStr[i] = ' ';
      break;
    }
  }
  atomInfo.residueId =  stoi(residueIdStr);


  // complete alternative residue determination
  // all conclude in residue name
  if(altTag != 'X'){
    atomInfo.residue.insert(atomInfo.residue.begin(), altTag);
  }



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

  if(atomInfo.element == 'H') return;
  // ignore H atoms
  inferChain(atomInfo.chainNum, atomInfo.residueId);  // make sure chain exists
  atomInfo.residueId += chains[atomInfo.chainNum].offset;
  if (lastResId != atomInfo.residueId){
    // a new residue now, checkout the last residue
    if(lastResId != -1) cleanAltRes();
    lastResId = atomInfo.residueId;
  }
  // push this one to altRes buffer
  int index;
  for(index = 0; index < altRes.size(); index++){
    // check whether the residue already exists in the buffer
    if(altRes[index].type == atomInfo.residue){
      break;
    }
  }
  if(index == altRes.size()){
    // doesn't exist, create a new residue
    Residue r(atomInfo.residue, atomInfo.residueId, atomInfo.chainNum);
    altRes.push_back(r);
  }
  if(atomInfo.residue.length() > 3) atomInfo.residue.erase(0, 1);
  Atom atom(atomInfo);
  // keep atom clean
  altRes[index].atoms.push_back(atom);
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
  // cout << "residue ID" << residueId << endl;
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
  cout << "Writing " << fileName << endl;
  fout.open(fileName.c_str());
  Chain & chain = chains[chainNum];
  vector<Residue> &residues = chain.residues;
  int lineNum = 1;
  bool baseIndexSet = false;
  int baseIndex = 0;
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
          fout << "  "; // for chainNum;
          if(!baseIndexSet){
             baseIndex = atom.residueId;
             baseIndexSet = true;
          }
          int oneBasedIndex = atom.residueId - baseIndex + 1;
          fout << setw(4) << oneBasedIndex << "    ";// change to 1 based
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
  fout << "TER\n";
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