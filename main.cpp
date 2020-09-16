#include <iostream>

#include "Brisk2.hpp"

using namespace std;


int main(int argc, char** argv) {
	if (argc < 2) {
		cout << "[fasta file]" << endl;
		cout << "0 DEFAULT, OUTPUT THE COUNT WITH NO CHECKING" << endl;
		cout << "1 PERFORMANCE MODE NO COUNT OUTPUT" << endl;
		cout << "2 DEBUG MODE COUNT AND CHECK WITH HASHTABLE" << endl;
		exit(0);
	}
	
	bool check;
	int mode(2);
	if (argc > 2) {
		mode = stoi(argv[2]);
	}

	if (mode > 1) {
		check = true;
		cout << "LETS CHECK THE RESULTS" << endl;
	}


	cout << "\n\n\nI count " << argv[1] << endl;
	cout << "Minimizer size:	" << minimizer_size << endl;
	cout << "Kmer size:	" << k << endl;
  
	auto start = std::chrono::system_clock::now();
	read_fasta_file(argv[1]);
	
	auto end = std::chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	cout << "Kmer counted elapsed time: " << elapsed_seconds.count() << "s\n";
	cout << endl;

	if (mode == 2) {
		cout<<menu.dump_counting()<<" errors"<<endl;;
	}
	menu.dump_stats();
}
