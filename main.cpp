#include <iostream>
#include <stdint.h>

#include "CLI11.hpp"
#include "Brisk2.hpp"

using namespace std;


int parse_args(int argc, char** argv, string & fasta, uint8_t & k, uint8_t & m,
								uint & mode) {
	CLI::App app{"Brisk library demonstrator - kmer counter"};

  auto file_opt = app.add_option("-f,--file", fasta, "Fasta file to count");
  file_opt->required();
  app.add_option("-k", k, "Kmer size");
  app.add_option("-m", m, "Minimizer size");
  app.add_option("--mode", mode, "Execution mode (0: output count, no checking | 1: performance mode, no output | 2: debug mode");

  CLI11_PARSE(app, argc, argv);
  return 0;
}


int main(int argc, char** argv) {
	string fasta = "";
	uint8_t k=31, m=11;
	uint mode = 0;

	parse_args(argc, argv, fasta, k, m, mode);
  cout << fasta << " " << (uint)k << " " << (uint)m << endl;

	// if (mode > 1) {
	// 	check = true;
	// 	cout << "LETS CHECK THE RESULTS" << endl;
	// }


	cout << "\n\n\nI count " << fasta << endl;
	cout << "Kmer size:	" << (uint)k << endl;
	cout << "Minimizer size:	" << (uint)m << endl;
  
	// auto start = std::chrono::system_clock::now();
	// read_fasta_file(argv[1]);
	
	// auto end = std::chrono::system_clock::now();
	// chrono::duration<double> elapsed_seconds = end - start;
	// cout << "Kmer counted elapsed time: " << elapsed_seconds.count() << "s\n";
	// cout << endl;

	// if (mode == 2) {
	// 	cout<<menu.dump_counting()<<" errors"<<endl;;
	// }
	// menu.dump_stats();

	return 0;
}
