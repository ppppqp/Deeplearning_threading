// This is a testing file
#include <stdlib.h>
#include <unistd.h>

#include <fstream>
#include <sstream>
#include <string>

#include "fileParser.cpp"
int main() {
  std::ifstream reference;
  reference.open("../../summary/AAA.seq");
  string line;
  double successCount = 0;
  double failureCount = 0;
  while (getline(reference, line)) {
    if (line[0] == '>') {  // start with file indicator
      std::stringstream ss(line);
      string indicator;
      int resNum;
      ss >> indicator >> resNum;
      if (indicator.size() > 6) continue;  // skip sub chain
      char chainNum = indicator[5];
      string proteinName = indicator.substr(1, 5);
      string pdbFileName = "../../PDB/" + proteinName + ".pdb";
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
      PDBParser parser(pdbFileName, true);

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
      } else {
        successCount++;
      }
    }
  }
  cout << "Passed: " << successCount << "\nFailed: "<< failureCount<< '\n';
}
