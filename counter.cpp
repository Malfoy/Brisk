#include <iostream>
#include <stdint.h>
#include <chrono>

#include "CLI11.hpp"
#include "zstr.hpp"
#include "brisk/Brisk.hpp"
#include "brisk/writer.hpp"


using namespace std;

// --- Useful functions to count kmers ---
void count_fasta(Brisk<uint8_t> & counter, string & filename, const uint threads);
void count_sequence(Brisk<uint8_t> & counter, string & sequence);
void verif_counts(Brisk<uint8_t> & counter);


int parse_args(int argc, char** argv, string & fasta, string & outfile, uint8_t & k, uint8_t & m, uint8_t & buckets,
								uint & mode, uint & threads) {
	CLI::App app{"Brisk library demonstrator - kmer counter"};

  auto file_opt = app.add_option("-f,--file", fasta, "Fasta file to count");
  file_opt->required();
  app.add_option("-k", k, "Kmer size");
  app.add_option("-m", m, "Minimizer size");
  app.add_option("-b", buckets, "Bucket order of magnitude. 4^b minimizer per bucket");
  app.add_option("-t", threads, "Thread number");
  app.add_option("-o", outfile, "Output file (kff format https://github.com/yoann-dufresne/kmer_file_format)");
  app.add_option("--mode", mode, "Execution mode (0: output count, no checking | 1: performance mode, no output | 2: debug mode");

  CLI11_PARSE(app, argc, argv);

  return 0;
}


string pretty_int(uint64_t n){
	string result;
	uint64_t order=1000000000000000000;
	bool started(false);
	while(order!=1){
		if(n/order>=1){
			string local( to_string(n/order));
			if(started){
				if(local.size()==2){
					result+='0';
				}
				if(local.size()==1){
					result+="00";
				}
			}
			result+=local+",";
			started=true;
			n%=order;
		}else if (started){
			result+="000,";
		}
		order/=1000;
	}
	string local( to_string(n));
	if(started){
		if(local.size()==2){
			result+='0';
		}
		if(local.size()==1){
			result+="00";
		}
	}
	result+=local;

	return result;
}

static robin_hood::unordered_map<kint, int16_t> verif;
static bool check;
static uint64_t number_kmer_count(0);


int main(int argc, char** argv) {
	string fasta = "";
	string outfile = "";
	uint8_t k=63, m=13, b=4;
	uint mode = 0;
	uint threads = 8;

	if (parse_args(argc, argv, fasta, outfile, k, m, b, mode, threads) != 0 or fasta == "")
		exit(0);

	Parameters params(k, m, b);
  cout << fasta << " " << (uint)k << " " << (uint)m << endl;

	if (mode > 1) {
		check = true;
		cout << "LETS CHECK THE RESULTS" << endl;
	}


	cout << "\n\n\nI count " << fasta << endl;
	cout << "Kmer size:	" << (uint)k << endl;
	cout << "Minimizer size:	" << (uint)m << endl;
  
	auto start = std::chrono::system_clock::now();
	Brisk<uint8_t> counter(params);
	count_fasta(counter, fasta, threads);
	
	auto end = std::chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	cout << "Kmer counted elapsed time: " << elapsed_seconds.count() << "s\n";
	cout << endl;


	if (check)
		verif_counts(counter);

	cout << "Global statistics:" << endl;
	uint64_t nb_buckets, nb_skmers, nb_kmers, nb_cursed, memory;
	counter.stats(nb_buckets, nb_skmers, nb_kmers, nb_cursed, memory);
	cout << pretty_int(nb_buckets) << " bucket used (/" << pretty_int(pow(4, counter.params.m_small)) << " possible)" << endl;
	cout << "nb superkmers: " << pretty_int(nb_skmers) << endl;
	cout << "nb kmers: " << pretty_int(nb_kmers) << endl;
	cout << "kmer / second: " << ((float)nb_kmers / elapsed_seconds.count()) << endl;
	cout << "average kmer / superkmer: " << ((float)nb_kmers / (float)nb_skmers) << endl;
	cout << "average superkmer / bucket: " << ((float)nb_skmers / (float)nb_buckets) << endl;
	cout << "Memory usage: " << (memory / 1024) << "Mo" << endl;
	cout << "bits / kmer: " << ((float)(memory * 1024 * 8) / (float)nb_kmers) << endl;
	cout << "nb cursed kmers: " << pretty_int(nb_cursed) << endl;
	cout<<"nb kmer considered: " <<pretty_int(number_kmer_count)<<endl;

	// --- Save Brisk index ---
	if (mode == 0 and outfile != "") {
		BriskWriter writer(outfile);
		writer.write(counter);
		writer.close();
	}

	return 0;
}


void verif_counts(Brisk<uint8_t> & counter) {
	cout << "--- Start counting verification ---" << endl;

	// kint mini_mask = (1 << (2 * counter.m)) - 1;
	kmer_full kmer(0,0, counter.params.m, false);
	// Count 
	while (counter.next(kmer)) {
		if (verif.count(kmer.kmer_s) == 0) {
			verif[kmer.kmer_s] = 0;
		}


		uint8_t * count = counter.get(kmer);
		if (count == NULL) {
			print_kmer(kmer.minimizer, counter.params.m); cout << endl;
			cout << (uint)kmer.minimizer << endl;
			cout << (uint*)count << endl;
			print_kmer(kmer.kmer_s, counter.params.k); cout << endl;
			cout << *(((uint*)&(kmer.kmer_s))+1) << " " << (uint)kmer.kmer_s << endl;
		}
		verif[kmer.kmer_s] -= *count;

		kmer.kmer_s = 0;
	}

	// Summary print
	uint errors = 0;
	for (auto & it : verif) {
		if (it.second != 0) {
			errors += 1;
			if (it.second > 0) {
				cout << "missing "; print_kmer(it.first, counter.params.k); cout << " " << (uint)it.second << endl;
			} else {
				cout << "too many "; print_kmer(it.first, counter.params.k); cout << " " << (uint)(-it.second) << endl;
				// cout << (uint)it.first << endl;
			}
		}
	}

	if (errors == 0){
		cout << "All counts are correct !" << endl;
	}else{
		cout<<errors << " errors" << endl;
	}

	cout << endl;
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
			case 'N':break;
			case 'n':break;
			default: str[i]='A';
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
void count_fasta(Brisk<uint8_t> & counter, string & filename, const uint threads) {
	// Test file existance
	struct stat exist_buffer;
  bool file_existance = (stat (filename.c_str(), &exist_buffer) == 0);
	if(not file_existance){
		cerr<<"Problem with file opening:	"<<filename<<endl;
		exit(1);
	}

	// Read file line by line
	cout << filename << " " << filename.length() << endl;
	zstr::ifstream in(filename);
	vector<string>  buffer;

	// uint idx = 0;
	#pragma omp parallel num_threads(threads)
	{
		while (in.good() or not buffer.empty()) {
			string line;
			#pragma omp critical
			{
				if(not buffer.empty()){
					line=buffer[buffer.size()-1];
					buffer.pop_back();
				}else{
					if (in.good()) {
						line = getLineFasta(&in);
						if(line.size()>1000000){
							buffer.push_back(line.substr(0,line.size()/4));
							buffer.push_back(line.substr(line.size()/4-counter.params.k+1,line.size()/4+counter.params.k-1));
							buffer.push_back(line.substr(line.size()/2-counter.params.k+1,line.size()/4+counter.params.k-1));
							buffer.push_back(line.substr(3*line.size()/4-counter.params.k+1));
						}
					} else
						line = "";
				}
			}

			if (line != "") {
				count_sequence(counter, line);
				// #pragma omp critical
				// {
				// 	if (idx % 10000 == 0)
				// 		cout << idx << endl;
				// 	idx++;
				// }
			}
		}
	}
}



void count_sequence(Brisk<uint8_t> & counter, string & sequence) {
	// Line too short
	if (sequence.size() < counter.params.k)
		return;

	vector<vector<kmer_full> > kmers_by_minimizer;
	vector<kmer_full> superkmer;

	// kint minimizer;
	SuperKmerEnumerator enumerator(sequence, counter.params.k, counter.params.m);

	// minimizer = enumerator.next(superkmer);
	enumerator.next(superkmer);
	while (superkmer.size() > 0) {
		// Add the values
		for (kmer_full & kmer : superkmer) {
			if (check) {
				#pragma omp critical
				{
					if (verif.count(kmer.kmer_s) == 0){
						verif[kmer.kmer_s] = 0;
					}
						
					verif[kmer.kmer_s] += 1;
					verif[kmer.kmer_s] = verif[kmer.kmer_s] % 256;
				}
			}
			if(kmer.suffix_size()>=0){
				#pragma omp atomic
				number_kmer_count++;
				counter.protect_data(kmer);
				uint8_t * data_pointer = counter.get(kmer);

				if (data_pointer == NULL) {
					data_pointer = counter.insert(kmer);
					// Init counter
					*data_pointer = (uint8_t)0;
				}
				// Increment counter
				*data_pointer += 1;
				counter.unprotect_data(kmer);
			}

			
		}

		// Next superkmer
		superkmer.clear();
		// minimizer = enumerator.next(superkmer);
		enumerator.next(superkmer);
	}

	return;
}
