// This is a testing file
#include <stdlib.h>
#include <unistd.h>

#include <fstream>
#include <ios>
#include <iostream>
#include <sstream>
#include <string>

#include "fileParser.cpp"

int main() {
  std::ifstream reference;
  std::ifstream pdbList;
  std::ifstream completePDB;
  double vm, rss;
  try {
    pdbList.open("../data/pdb_list");
    string line;
    getline(pdbList, line);// table header;
    while(getline(pdbList, line)){
        string proteinName = line.substr(0, 4);
        char chainNum = line[4];
        string chainNumStr = line.substr(4,1);
        string completePDBFile = "../data/pdb/"+proteinName+".pdb";
        if (access(completePDBFile.c_str(), F_OK) != 0) {
          // no avaiable pdb file
          continue;
        }
        // completePDB.open(completePDBFile);
        // if(!completePDB.is_open()){
        //     cout << "fail to open file" + completePDBFile << endl;
        // }

        // std::transform(proteinName.begin(), proteinName.end(),
        //         proteinName.begin(), ::tolower);
        // string referenceFile = "../../PDB/"+ proteinName + chainNumStr+".pdb";
        // reference.open(referenceFile);
        // if(!reference.is_open()){
        //     cout << "fail to open file" + referenceFile << endl;
        //     continue;
        // }
        // cout << "successfully open files!" << endl;
        PDBParser parser(completePDBFile);
        parser.parse();
        parser.output2PDB(chainNum, "../data/clean/");
    }
    /*
    double successCount = 0;
    double failureCount = 0;
    int lineCount = 0;
    while (getline(reference, line)) {
      lineCount++;
      if (line[0] == '>') {  // start with file indicator
        std::stringstream ss(line);
        string indicator;
        int resNum;
        ss >> indicator >> resNum;
        if (indicator.size() > 6) continue;  // skip sub chain
        char chainNum = indicator[5];
        string proteinName = indicator.substr(1, 4);
        std::transform(proteinName.begin(), proteinName.end(),
                       proteinName.begin(), ::toupper);
        // cout << "protein name" << proteinName << endl;
        string pdbFileName = "../data/pdb/" + proteinName + ".pdb";
        if (access(pdbFileName.c_str(), F_OK) != 0) {
          // no avaiable pdb file
          continue;
        }
        // start parsing correct answer
        string correct;
        int lineNum = ceil(static_cast<double>(resNum) / 70);
        for (int i = 0; i < lineNum; i++) {
          getline(reference, line);
          correct += line;
        }
        // using pdbParser to generat fasta file
        PDBParser parser(pdbFileName);

        parser.parse();

        parser.output2Fasta();
        // start parsing output file
        std::ifstream outputFile;

        outputFile.open("pdb2fasta");
        string output;

        while (getline(outputFile, line)) {
          if (line[0] == '>') {
            std::stringstream outputSS(line);
            string outputIndicator;
            int outputResNum;

            outputSS >> outputIndicator >> outputResNum;
            if (outputIndicator == indicator) {
              int outputLineNum = ceil(static_cast<double>(outputResNum) / 70);
              for (int i = 0; i < outputLineNum; i++) {
                getline(outputFile, line);
                output += line;
              }
            }
          }
        }
        if (output != correct) {
          failureCount++;
          cout << setw(30) << "Test protein " << proteinName << " chain "
               << chainNum << " test FAILED!\n";
          cout << "output is  " << output << '\n';
          cout << "correct is " << correct << '\n';
          for(int i = 0; i < output.size(); i++) if(output[i]!=correct[i]){
            cout << "first index: " << i << endl;
            break;
          }
        } else {
          successCount++;
          cout << "Test protein " << proteinName << " chain " << chainNum
               << " test SUCCESS! "
               << "Passed: " << successCount << " Failed: " << failureCount
               << '\n';
          ;
        }
      }
    }
    */
  }
  catch(char const * error){
    cout << error << endl;
  }
}