output: main.cpp fileParser.cpp Atom.cpp Chain.cpp Residue.cpp
	g++ -std=c++11 main.cpp -o output

test: test_fasta.cpp fileParser.cpp Atom.cpp Chain.cpp Residue.cpp
	g++ -std=c++11 test_fasta.cpp -o test

test_linux: test_fasta_linux.cpp fileParser.cpp Atom.cpp Chain.cpp Residue.cpp
	g++ -std=c++11 test_fasta_linux.cpp -o test_linux

test_data: test_fasta_data.cpp fileParser.cpp Atom.cpp Chain.cpp Residue.cpp
	g++ -std=c++11 test_fasta_data.cpp -o test_data

test_pdb: test_writePDB.cpp fileParser.cpp Atom.cpp Chain.cpp Residue.cpp
	g++ -std=c++11 test_writePDB.cpp -o test_pdb

clean: 
	rm *.o output test_linux test test_data


# g++ -o  test -std=c++11 test.cpp