output: main.cpp
	g++ -std=c++11 main.cpp -o output

test: test_fasta.cpp
	g++ -std=c++11 test_fasta.cpp -o test

test_linux: test_fasta_linux.cpp
	g++ -std=c++11 test_fasta_linux.cpp -o test_linux

test_data: test_fasta_data.cpp
	g++ -std=c++11 test_fasta_data.cpp -o test_data

clean: 
	rm *.o output test_linux test test_data


# g++ -o  test -std=c++11 test.cpp