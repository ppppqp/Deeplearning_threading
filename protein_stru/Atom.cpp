#ifndef Atom_cpp
#define Atom_cpp
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using std::cout;
using std::endl;
using std::setw;
using std::string;

const int LINEWIDTH = 23;

class AtomInfo {
 public:
  string atomType;
  string residue;
  double x, y, z, occupancy, beta;
  int atomId;
  int residueId;
  char element;
  char chainNum;
};
class Atom {
 public:
  string atomType;
  string residue;
  double x, y, z, occupancy, beta;
  int atomId;
  int residueId;
  char element;
  char chainNum;
  bool valid;
  /**
   * Atom Class Default Ctor
   **/
  Atom(bool isValid)
      : atomId(0),
        atomType(""),
        residue(""),
        chainNum('X'),
        residueId(0),
        x(0),
        y(0),
        z(0),
        occupancy(0),
        beta(0),
        element('X'),
        valid(isValid){};
  Atom()
      : atomId(0),
        atomType(""),
        residue(""),
        chainNum('X'),
        residueId(0),
        x(0),
        y(0),
        z(0),
        occupancy(0),
        beta(0),
        element('X'),
        valid(true){};
  /**
   * Atom Class Constructor from argument
   **/
  Atom(int _atomId, string _atomType, string _residue, int _chainNum,
       int _residueId, double _x, double _y, double _z, double _occupancy,
       double _beta, char _element)
      : atomId(_atomId),
        atomType(_atomType),
        residue(_residue),
        chainNum(_chainNum),
        residueId(_residueId),
        x(_x),
        y(_y),
        z(_z),
        occupancy(_occupancy),
        beta(_beta),
        element(_element),
        valid(true){

        };
  /**
   * Atom Class Ctor from AtomInfo Object
   **/
  Atom(AtomInfo atomInfo)
      : atomId(atomInfo.atomId),
        atomType(atomInfo.atomType),
        residue(atomInfo.residue),
        chainNum(atomInfo.chainNum),
        residueId(atomInfo.residueId),
        x(atomInfo.x),
        y(atomInfo.y),
        z(atomInfo.z),
        occupancy(atomInfo.occupancy),
        beta(atomInfo.beta),
        element(atomInfo.element),
        valid(true){};
  char getAtomId() { return atomId; }
  string getAtomType() { return atomType; }
  char getChainNum() { return chainNum; }
  string getResidue() { return residue; }
  double getX() { return x; }
  double getY() { return y; }
  double getZ() { return z; }
  double getOccupancy() { return occupancy; };
  double getBeta() { return beta; };
  char getElement() { return element; };
  friend std::stringstream &operator>>(std::stringstream &ss, Atom &A);
  void print() {
    cout << "###" << setw(LINEWIDTH) << "\t\t---ATOM INFO START---" << endl;
    cout << "\t a" << setw(LINEWIDTH) << "\t\tatomId: " << atomId << endl;
    cout << "\t a" << setw(LINEWIDTH) << "\t\tatomType: " << atomType << endl;
    cout << "\t a" << setw(LINEWIDTH) << "\t\tresidue: " << residue << endl;
    cout << "\t a" << setw(LINEWIDTH) << "\t\tchainNum: " << chainNum << endl;
    cout << "\t a" << setw(LINEWIDTH) << "\t\tresidueId: " << residueId << endl;
    cout << "\t a" << setw(LINEWIDTH) << "\t\tx: " << x << endl;
    cout << "\t a" << setw(LINEWIDTH) << "\t\ty: " << y << endl;
    cout << "\t a" << setw(LINEWIDTH) << "\t\tz: " << z << endl;
    cout << "\t a" << setw(LINEWIDTH) << "\t\toccupancy: " << occupancy << endl;
    cout << "\t a" << setw(LINEWIDTH) << "\t\tbeta: " << beta << endl;
    cout << "\t a" << setw(LINEWIDTH) << "\t\telement: " << element << endl;
    cout << "###" << setw(LINEWIDTH) << "\t\t---ATOM INFO END---" << endl;
  }
};
#endif