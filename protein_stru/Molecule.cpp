#ifndef Molecule_cpp
#define Molecule_cpp
#include <map>
#include <string>
#include <vector>
using std::string;
using std::vector;
#include "Atom.cpp"
#include "Chain.cpp"
using std::map;
using std::setw;
class Molecule {
 public:
  vector<char> chainIdList;
  int id;
  string name;
  string synonym;
  string engineered;
  string scientific;
  string common;
  int taxId;
  vector<string> gene;
  void print() {
    cout << setw(30) << "---MOLECULE INFO START---" << endl;
    cout << setw(30) << "id: " << id << endl;
    cout << setw(30) << "name: " << name << endl;
    cout << setw(30) << "---MOLECULE INFO END---" << endl;
  }
};
#endif