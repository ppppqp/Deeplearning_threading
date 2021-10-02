output: Atom.o fileParser.o main.o
	g++ -std=c++11 Atom.o  main.o -o output

test: fileParser.o test_fasta.o
	g++ -std=c++11 test_fasta.o -o test

test_linux: fileParser.o test_fasta_linux.o
	g++ -std=c++11 test_fasta_linux.o -o test_linux

readChain.o: readChain.cpp
	g++ -c  readChain.cpp

test_fasta.o: test_fasta.cpp
	g++ -c test_fasta.cpp

test_fasta_linux.o: test_fasta_linux.cpp
	g++ -c test_fasta_linux.cpp

Residue.o: Residue.cpp
	g++ -c  Residue.cpp

Atom.o: Atom.cpp
	g++ -c  Atom.cpp

fileParser.o: fileParser.cpp
	g++ -c fileParser.cpp

Chain.o: Chain.cpp
	g++ -c  Chain.cpp

main.o: main.cpp
	g++ -c  main.cpp

clean: 
	rm *.o output