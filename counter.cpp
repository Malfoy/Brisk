#include <iostream>
#include <stdint.h>
#include <chrono>
#include <unordered_map>

#include "CLI11.hpp"
#include "zstr.hpp"
#include "brisk/buckets.hpp"
#include "brisk/Brisk.hpp"
#include "brisk/Kmers.hpp"
#include "brisk/writer.hpp"



using namespace std;



// --- Useful functions to count kmers ---
void count_fasta(Brisk<uint8_t> & counter, string & filename, const uint threads);
void count_sequence(Brisk<uint8_t> & counter, string & sequence);
void verif_counts(Brisk<uint8_t> & counter);



int parse_args(int argc, char** argv, string & fasta, string & outfile, uint8_t & k, uint & mode, uint & threads) {
	CLI::App app{"Brisk library demonstrator - kmer counter"};

  auto file_opt = app.add_option("-f,--file", fasta, "Fasta file to count")->required();
  file_opt->check(CLI::ExistingFile);
  app.add_option("-k", k, "Kmer size")->check(CLI::Range(5,63));
  app.add_option("-t", threads, "Number of threads to use");
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



// static robin_hood::unordered_map<kint, int16_t> verif;
// static tsl::sparse_map<kint, int16_t> verif;
// typename ankerl::unordered_dense::map<kint, int16_t> verif;
static unordered_map<kint, int16_t> verif;
static bool check;



int main(int argc, char** argv) {
	string fasta = "";
	string outfile = "";
	uint8_t k=63, m=11, b=6;
	uint mode = 0;
	uint threads = 8;

	if (parse_args(argc, argv, fasta, outfile, k, mode, threads) != 0 or fasta == "")
		exit(0);

	Parameters params(k, m, b);

	if (mode > 1) {
		check = true;
		cout << "LETS CHECK THE RESULTS" << endl;
	}


	cout << "I'm counting " << fasta << endl;
	cout << "Kmer size:	" << (uint)k << endl;
	cout << "Minimizer size:	" << (uint)m << endl;
	cout << "Bucket size:     " << (uint)b << endl;
  
	auto start = std::chrono::system_clock::now();
	Brisk<uint8_t> counter(params);	
	count_fasta(counter, fasta, threads);
	
	auto end = std::chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
	cout << "Kmer counted elapsed time: " << elapsed_seconds.count() << "s\n";
	cout << endl;

	if (check)
		verif_counts(counter);

	uint64_t nb_buckets, nb_skmers, nb_kmers, memory, largest_bucket;
	counter.stats(nb_buckets, nb_skmers, nb_kmers, memory, largest_bucket);
	cout << pretty_int(nb_buckets) << " bucket used (/" << pretty_int(pow(4, counter.params.b)) << " possible)" << endl;
	cout << "nb superkmers: " << pretty_int(nb_skmers) << endl;
	cout << "nb kmers: " << pretty_int(nb_kmers) << endl;
	cout << "kmer / second: " << pretty_int((float)nb_kmers / elapsed_seconds.count()) << endl;
	cout << "average kmer / superkmer: " << ((float)nb_kmers / (float)nb_skmers) << endl;
	cout << "average superkmer / bucket: " << ((float)nb_skmers / (float)nb_buckets) << endl;
	cout << "Largest bucket :	"<<pretty_int(largest_bucket) <<endl;
	cout << "Memory usage: " << (memory / 1024) << "Mo" << endl;
	cout << "bits / kmer: " << ((float)(memory * 1024 * 8) / (float)nb_kmers) << endl;

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
	kmer_full kmer(0,0, counter.params.m, counter.params.dede);
	// Count 
	while (counter.next(kmer)) {
		uint8_t * count = counter.get(kmer);
		if (count == NULL) {
			cout<<"No data linked ";
			print_kmer(kmer.kmer_s, counter.params.k); cout << endl;
		}else{
			verif[kmer.kmer_s] -= *count;
		}

		kmer.kmer_s = 0;
	}

	// Summary print
	uint errors = 0;
	for (auto & it : verif) {
		if (it.second != 0) {
			errors += 1;
			if (it.second > 0) {
				cout << "missing " << (uint)it.second << " "; print_kmer(it.first, counter.params.k);
			} else {
				cout << "too many " << (uint)(-it.second) << " "; print_kmer(it.first, counter.params.k);
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



void clean_dna(string& str,string& previous){
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
			default: 
			previous=str.substr(i+1);
			str=str.substr(0,i);
		}
	}
	transform(str.begin(), str.end(), str.begin(), ::toupper);
}



string getLineFasta(zstr::ifstream* in,string& previous) {
	string line, result;
	if(previous.empty()){
		getline(*in, line);
		char c = static_cast<char>(in->peek());
		while (c != '>' and c != EOF) {
			getline(*in, line);
			result += line;
			c = static_cast<char>(in->peek());
		}
	}else{
		result=previous;
		previous="";
	}
	
	clean_dna(result,previous);
	return result;
}



/** Counter function.
  * Read a complete fasta file line by line and store the counts into the Brisk datastructure.
  */
void count_fasta(Brisk<uint8_t> & counter, string & filename, const uint threads) {
	kmer_full init;
	// Test file existance
	ifstream fs;
	fs.open(filename);
  bool file_existance = fs ? true : false;
	if(not file_existance){
		cerr<<"Problem with file opening:	"<<filename<<endl;
		exit(1);
	}
  fs.close();

	// Read file line by line
	// cout << filename << " " << filename.length() << endl;
	zstr::ifstream in(filename);
	// omp_set_nested(2);
	//#pragma omp parallel
	{
		string line,prev;
		while (in.good() or prev.size()!=0) {
			
			#pragma omp critical
			{
				line = getLineFasta(&in,prev);
			}

			if (line != "") {
				count_sequence(counter, line);
			}
		}
	}
}



void count_sequence(Brisk<uint8_t> & counter, string & sequence) {
	// Line too short
	if (sequence.size() < counter.params.k){
		return;
	}
	omp_lock_t local_mutex;
	omp_init_lock(&local_mutex);
	SuperKmerEnumerator enumerator(sequence, counter.params.k, counter.params.m,counter.params.dede);
	// #pragma omp parallel
	{
		vector<kmer_full> superkmer;
		vector<bool> newly_inserted;
		vector<uint8_t*> vec;

		omp_set_lock(&local_mutex);
		kint minimizer = enumerator.next(superkmer);
		omp_unset_lock(&local_mutex);
		while (superkmer.size() > 0) {
			kmer_full local;
			local.copy(superkmer[0]);
			counter.protect_data(local);
			// Add the values
			if (check) {
				for (kmer_full & kmer : superkmer) {
					#pragma omp critical
					{
						if (verif.count(kmer.kmer_s) == 0){
							verif[kmer.kmer_s] = 0;
						}
						verif[kmer.kmer_s] += 1;
						verif[kmer.kmer_s] = verif[kmer.kmer_s] % 256;
					}
				}
			}
			newly_inserted.clear();

			vec=(counter.insert_superkmer(superkmer,newly_inserted));
			// cout << "newly_inserted size: " << newly_inserted.size() << endl;
			// cout << "vec size: " << vec.size() << endl;
			for(uint i(0); i < vec.size();++i){
				// cout << "vec[" << i << "]: " << &vec[i] << endl;
				uint8_t * data_pointer(vec[i]);
				// cout << "newly_inserted[" << i << "]: " << newly_inserted[i] << endl;
				if ((data_pointer)==NULL){
					cout << "PROBLEME" << endl; cin.get();
				}
				if(newly_inserted[i]){
					// cout << "1" << endl;
					(*data_pointer)=1;
				}else{
					// cout << "2" << endl;
					(*data_pointer)++;
				}
			}
			// cout << "3"  << endl;
			counter.unprotect_data(local);
			// Next superkmer
			superkmer.clear();
			// cout << "largest bucket: " << counter.menu->largest_bucket << endl;
			if (counter.menu->largest_bucket >= 65536) {
				cout << "Starting reallocation" << endl;
				counter.reallocate();
				cout << "Finished reallocation" << endl;
				enumerator.update(counter.params.m,counter.params.dede);
				string lanadine=sequence.substr(enumerator.seq_idx);
				enumerator.seq=lanadine;
				enumerator.seq_idx=0;
				// cin.get();	
			}

			omp_set_lock(&local_mutex);
			minimizer = enumerator.next(superkmer);
			omp_unset_lock(&local_mutex);

			if(minimizer==0){
				break;
			}
		}
	}
}
