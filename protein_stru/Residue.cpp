#ifndef Residue_cpp
#define Residue_cpp
#include <algorithm>
#include <vector>

#include "Atom.cpp"

struct compareAtomType {
  string atomType;
  compareAtomType(string const &i) : atomType(i) {}
  bool operator()(Atom const &i) { return (i.atomType == atomType); }
};
using std::vector;
class Residue {
 public:
  Residue();
  Residue(string _type, int _residueId, char _chainId);
  Residue(bool _valid);
  vector<Atom> atoms;  // atoms the residue contains
  int residueId;       // residue id
  char chainId;        // chain it belongs
  string type;
  bool valid;
  void print();
  Atom getAtom(string atomType);
};
Atom Residue::getAtom(string atomType) {
  for (int i = 0; i < atoms.size(); i++) {
    if (atoms[i].atomType == atomType) {
      Atom atom(atoms[i]);
      return atom;
    }
  }
  return Atom(false);
};
Residue::Residue() {}
Residue::Residue(bool _valid) : valid(_valid) {}
Residue::Residue(string _type, int _residueId, char _chainId)
    : type(_type),
      residueId(_residueId),
      chainId(_chainId),
      valid(true),
      atoms(vector<Atom>(0)) {}
void Residue::print() {
  cout << "###" << setw(20) << "\t--------RESIDUE INFO START--------" << endl;
  cout << 'r' << setw(20) << "\tresidueId: " << residueId << endl;
  cout << 'r' << setw(20) << "\ttype: " << type << endl;
  for (int i = 0; i < atoms.size(); i++) {
    // atoms[i].print();
  }
  cout << "###" << setw(20) << "\t--------RESIDUE INFO END--------" << endl;
}
#endif