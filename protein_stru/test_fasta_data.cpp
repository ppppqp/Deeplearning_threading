// This is a testing file
#include <stdlib.h>
#include <unistd.h>

#include <fstream>
#include <ios>
#include <iostream>
#include <sstream>
#include <string>

#include "fileParser.cpp"

void process_mem_usage(double& vm_usage, double& resident_set) {
  using std::ifstream;
  using std::ios_base;
  using std::string;

  vm_usage = 0.0;
  resident_set = 0.0;

  // 'file' stat seems to give the most reliable results
  //
  ifstream stat_stream("/proc/self/stat", ios_base::in);

  // dummy vars for leading entries in stat that we don't care about
  //
  string pid, comm, state, ppid, pgrp, session, tty_nr;
  string tpgid, flags, minflt, cminflt, majflt, cmajflt;
  string utime, stime, cutime, cstime, priority, nice;
  string O, itrealvalue, starttime;

  // the two fields we want
  //
  unsigned long vsize;
  long rss;

  stat_stream >> pid >> comm >> state >> ppid >> pgrp >> session >> tty_nr >>
      tpgid >> flags >> minflt >> cminflt >> majflt >> cmajflt >> utime >>
      stime >> cutime >> cstime >> priority >> nice >> O >> itrealvalue >>
      starttime >> vsize >> rss;  // don't care about the rest

  stat_stream.close();

  long page_size_kb = sysconf(_SC_PAGE_SIZE) /
                      1024;  // in case x86-64 is configured to use 2MB pages
  vm_usage = vsize / 1024.0;
  resident_set = rss * page_size_kb;
}

int main() {
  std::ifstream reference;
  double vm, rss;
  try {
    reference.open("../../summary/AAA.seq");
    // reference.open("./test.seq");
    string line;
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
  }
  catch(char const * error){
    cout << error << endl;
  }
}