#include <iostream>
#include <string>
#include <sstream>
#include <cstdlib>
#include "utils.hpp"
using namespace std;

int main(){
    // string line = "ATOM      1  N   MET A   1     157.673 191.424 197.711  1.00 31.81           N ";
    string line = "ATOM      1  N   MET     1     122.423 112.215 127.730";
    string atomType = line.substr(13,4);
    string residue = line.substr(17,4);
    char  chainNum = line[21];
    int residueId = stoi(line.substr(25,6).c_str());
    double x = stod(line.substr(31,7));
    double y = stod(line.substr(39, 7));
    double z = stod(line.substr(47,7));
    cout << "atomType" << atomType << endl;
    cout << "residue" << residue << endl;
    cout << "chain num" << chainNum << endl;
    cout << "residueId" << residueId << endl;
    cout << "x" << x << endl;
    cout << "y" << y << endl;
    cout << "z" << z << endl;
    double occupancy = -1;
    double beta = -1;
    char element = 'X';
    string eleStr = "";
    if(line.length() > 60) occupancy = stod(line.substr(56,4));
    if(line.length() > 66) beta = stod(line.substr(61,5));
    if(line.length() > 78 && line.substr(76,2)!="  ") {
        eleStr = line.substr(76,2);
        trim(eleStr);
        if(eleStr.length()!=0) element = eleStr[0];
    }
    else {
        if(!(atomType[0]>='0' && atomType[0] <= '9'))
            element = atomType[0];
        else element = atomType[1];
    }

    cout << "element" << element << endl;
    return 0;
}