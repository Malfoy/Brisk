#include <iostream>
#include <stdint.h>
#include <chrono>

#include "CLI11.hpp"
#include "zstr.hpp"
#include "Brisk2.hpp"
#include "Kmers.hpp"

using namespace std;

// --- Useful functions to count kmers ---
void count_fasta(Brisk<uint8_t> & counter, string & filename);
void count_sequence(Brisk<uint8_t> & counter, string & sequence);


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
	uint8_t k=63, m=9;
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
  
	auto start = std::chrono::system_clock::now();
	Brisk<uint8_t> counter(k, m);
	count_fasta(counter, fasta);
	
	auto end = std::chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	cout << "Kmer counted elapsed time: " << elapsed_seconds.count() << "s\n";
	cout << endl;

	// if (mode == 2) {
	// 	cout<<menu.dump_counting()<<" errors"<<endl;;
	// }
	// menu.dump_stats();

	return 0;
}


void clean_dna(string& str){
	for(uint i(0); i< str.size(); ++i){
		switch(str[i]){
			case 'a':break;
			case 'A':break;
			case 'c':break;
			case 'C':break;
			case 'g':break;
			case 'G':break;
			case 't':break;
			case 'T':break;
			// case 'N':break;
			// case 'n':break;
			default:  str[i]='A';break;;
		}
	}
	transform(str.begin(), str.end(), str.begin(), ::toupper);
}

string getLineFasta(zstr::ifstream* in) {
	string line, result;
	getline(*in, line);
	char c = static_cast<char>(in->peek());
	while (c != '>' and c != EOF) {
		getline(*in, line);
		result += line;
		c = static_cast<char>(in->peek());
	}
	clean_dna(result);
	return result;
}

/** Counter function.
  * Read a complete fasta file line by line and store the counts into the Brisk datastructure.
  */
void count_fasta(Brisk<uint8_t> & counter, string & filename) {
	// Test file existance
	struct stat exist_buffer;
  bool file_existance = (stat (filename.c_str(), &exist_buffer) == 0);
	if(not file_existance){
		cerr<<"Problem with file opening:	"<<filename<<endl;
		exit(1);
	}

	// Read file line by line
	uint nb_core = 1;
	zstr::ifstream in(filename);
	vector<string>  buffer;
	uint line_count = 0;

	#pragma omp parallel num_threads(nb_core)
	{
		string line;
		while (in.good() or not buffer.empty()) {
			#pragma omp critical(input)
			{
				if(not buffer.empty()){
					line=buffer[buffer.size()-1];
					buffer.pop_back();
				}else{
					line = getLineFasta(&in);
					if(line.size()>100000000000){
						buffer.push_back(line.substr(0,line.size()/4));
						buffer.push_back(line.substr(line.size()/4-counter.k+1,line.size()/4+counter.k-1));
						buffer.push_back(line.substr(line.size()/2-counter.k+1,line.size()/4+counter.k-1));
						line=line.substr(3*line.size()/4-counter.k+1);
					}
				}

			}
			count_sequence(counter, line);
			line_count++;
		}
	}
}


void count_sequence(Brisk<uint8_t> & counter, string & sequence) {
	// Line too short
	if (sequence.size() < counter.k)
		return;

	vector<vector<kmer_full> > kmers_by_minimizer;
	vector<kmer_full> superkmer;

	string_to_kmers_by_minimizer(sequence, superkmer, counter.k, counter.m);
	while (superkmer.size() > 0) {
		// TODO: DO something here
		// cout << superkmer.size() << endl;
		// for (kmer_full & kmer : superkmer) {
		// 	kmer.print(counter.k, counter.m);
		// }

		superkmer.clear();
		string_to_kmers_by_minimizer(sequence, superkmer, counter.k, counter.m);
	}

	return;

	for (auto kmers: kmers_by_minimizer) {
		uint8_t * data_array = counter.insert(kmers);
	}
}
