#ifndef Chain_cpp
#define Chain_cpp
#include <vector>

#include "Residue.cpp"
#include "Vector.cpp"
using std::vector;
class Chain {
 public:
  int id;
  string name;
  string synonym;
  string engineered;
  string scientific;
  string common;
  int validResNum;
  int taxId;
  char chainNum;             // chain number
  vector<Residue> residues;  // the residue it contains
  Chain();
  Chain(Chain& c);
  Chain(char _chainNum);
  void L1Depth(string mode, string pattern, double localR);
  Chain getCAChain();
  void print();
};
Chain::Chain() {}
Chain::Chain(char _chainNum)
    : chainNum(_chainNum), validResNum(0), residues(vector<Residue>()){};
Chain::Chain(Chain& c)
    : chainNum(c.chainNum),
      id(c.id),
      name(c.name),
      synonym(c.synonym),
      engineered(c.engineered),
      scientific(c.scientific),
      common(c.common),
      taxId(c.taxId),
      residues(c.residues){};
Chain Chain::getCAChain() {
  Chain newChain(chainNum);
  for (int i = 0; i < residues.size(); i++) {
    Residue& res = residues[i];
    Atom CA = res.getAtom("CA");
    if (CA.valid) {
      Residue newRes(res.type, res.residueId, res.chainId);
      newRes.atoms.clear();
      newRes.atoms.push_back(CA);
      newChain.residues.push_back(newRes);
    }
  }
  return newChain;
}

void Chain::print() {
  cout << setw(30) << "--------CHAIN INFO START--------" << endl;
  cout << setw(30) << "chainNum: " << chainNum << endl;
  cout << setw(30) << "id: " << id << endl;
  cout << setw(30) << "name: " << name << endl;
  for (int i = 0; i < residues.size(); i++) {
    cout << setw(5) << residues[i].residueId << ':';
    cout << residues[i].type << ' ';
    if ((i + 1) % 13 == 0) cout << '\n';
  }
  cout << '\n';
  for (int i = 0; i < residues.size(); i++) {
    //   residuesi].print();
  }
  cout << setw(30) << "--------CHAIN INFO END--------" << endl;
}
void Chain::L1Depth(string mode, string pattern, double localR = 99999) {
  if (mode == "") mode = "atom-atom";
  string firstMode = mode.substr(0, mode.find('-'));
  string lastMode = mode.substr(mode.find('-') + 1);
  vector<Residue> resFirst;
  vector<Residue> resLast;
  cout << "firstmode:" << firstMode << " lastmode: " << lastMode << endl;
  if (firstMode == "atom") {
    resFirst = residues;
  } else if (firstMode == "residue") {
    resFirst = getCAChain().residues;
  } else {
    cout << "Mode first argument is wrong\n";
    exit(1);
  }
  if (lastMode == "atom") {
    resLast = residues;
  } else if (lastMode == "residue") {
    resLast = getCAChain().residues;
  } else {
    cout << "Mode last argument is wrong\n";
    exit(1);
  }
  for (int i = 0; i < resFirst.size(); i++) {
    if (!resFirst[i].valid) continue;
    for (int j = 0; j < resFirst[i].atoms.size(); j++) {
      Atom atomFirst = resFirst[i].atoms[j];
      Vector firstV(atomFirst);
      Vector subV(0, 0, 0);
      int N = 0;
      for (int m = 0; m < resLast.size(); m++) {
        if (!resLast[m].valid) continue;
        for (int n = 0; n < resLast[m].atoms.size(); n++) {
          //   cout << "i=" << i << "j=" << j << "m=" << m << "n=" << n << endl;
          Atom atomLast = resLast[m].atoms[n];
          Vector lastV(atomLast);
          Vector tempV = lastV.subtract(firstV);
          tempV.normalize();
          if (pattern == "global") {
            subV = subV.plus(tempV);
            N++;
          } else if (pattern == "local") {
            if (tempV.norm() <= localR) {
              subV = subV.plus(tempV);
              N++;
            } else {
              cout << "Pattern argument is wrong\n";
              exit(1);
            }
          }
        }
      }
      // subV.print();
      double result = 1 - subV.norm() / N;
      // cout << "N = " << N << '\n';
      cout.precision(16);
      cout << std::fixed << result << '\n';
    }
  }
}
#endif