#include <cmath>
#include <iostream>

#include "Atom.cpp"
using std::fixed;
using std::vector;
class Vector {
 public:
  double x;
  double y;
  double z;
  Vector(double _x, double _y, double _z);
  Vector(Atom a);
  double norm();
  Vector plus(Vector v);
  Vector subtract(Vector v);
  Vector dot(Vector v);
  void print();
  void normalize();
};

Vector::Vector(double _x, double _y, double _z) : x(_x), y(_y), z(_z) {}
Vector::Vector(Atom a) : x(a.x), y(a.y), z(a.z){};
double Vector::norm() { return sqrt(x * x + y * y + z * z); }
Vector Vector::plus(Vector v) {
  Vector vc(x + v.x, y + v.y, z + v.z);
  return vc;
}
Vector Vector::subtract(Vector v) {
  Vector vc(x - v.x, y - v.y, z - v.z);
  return vc;
}
Vector Vector::dot(Vector v) {
  Vector vc(x * v.x, y * v.y, z * v.z);
  return vc;
}
void Vector::print() {
  cout.precision(15);
  cout << fixed << "( " << x << ", " << y << ", " << z << " )\n";
}
void Vector::normalize() {
  double norm = this->norm();
  if (norm == 0)
    return;
  else {
    x /= norm;
    y /= norm;
    z /= norm;
  }
}