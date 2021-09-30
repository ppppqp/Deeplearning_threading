output: Atom.o fileParser.o main.o
	g++ -std=c++11 Atom.o  main.o   -o output

test: fileParsor.o test.o
	g++ fileParsor.o test.o -o test

test.o: test.cpp
	g++ -c test.cpp

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