#include <iostream>
#include <stdint.h>
#include <chrono>
#include <unordered_map>

#include "CLI11.hpp"
#include "zstr.hpp"
#include "Brisk.hpp"
#include "writer.hpp"


using namespace std;

// --- Useful functions to count kmers ---
void add_color_from_fasta(Brisk<uint64_t> & colors, string & filename, const uint threads, uint64_t color);
void add_color_from_sequence(Brisk<uint64_t> & colors, string & sequence, uint64_t color);
void verif_colors(Brisk<uint64_t> & colors);


string getLineFasta(zstr::ifstream* in);


int parse_args(int argc, char** argv, string & fasta_list, string & outfile, uint8_t & k, uint8_t & m, uint8_t & buckets,
								uint & mode, uint & threads) {
	CLI::App app{"Brisk library demonstrator - kmer counter"};

  auto file_opt = app.add_option("-f,--file-list", fasta_list, "Text file containing the path of sample fasta files. The fasta paths should be one per line.");
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

static unordered_map<kint, int16_t> verif;
static bool check;
static uint64_t number_kmer_count(0);
uint64_t low_complexity_kmer(0);


int main(int argc, char** argv) {
	string fasta_list = "";
	string outfile = "";
	uint8_t k=63, m=13, b=4;
	uint mode = 0;
	uint threads = 8;

	if (parse_args(argc, argv, fasta_list, outfile, k, m, b, mode, threads) != 0 or fasta_list == "")
		exit(0);

	Parameters params(k, m, b);
  cout << fasta_list << " " << (uint)k << " " << (uint)m << endl;

	if (mode > 1) {
		check = true;
		cout << "LETS CHECK THE RESULTS" << endl;
	}

	// Open the file containing the paths to each fasta of the collection
	ifstream fastas(fasta_list);
	if (!fastas.is_open()) {
		cerr << "Could not open the file - '" << fasta_list << "'" << endl;
		return EXIT_FAILURE;
	}

	auto global_start = std::chrono::system_clock::now();
	cout << "Kmer size:	" << (uint)k << endl;
	cout << "Minimizer size:	" << (uint)m << endl;

	Brisk<uint64_t> colors(params);
	uint64_t color_id = 1;
	string fasta;
	while (getline(fastas, fasta))
	{
		if (fasta.length() == 0)
			continue;

		cout << "\n\n\nI color " << fasta << endl;
	  
		auto start = std::chrono::system_clock::now();
		add_color_from_fasta(colors, fasta, threads, color_id);
		
		auto end = std::chrono::system_clock::now();
		chrono::duration<double> elapsed_seconds = end - start;
		cout << "Color elapsed time: " << elapsed_seconds.count() << "s\n";
		cout << endl;

		color_id <<= 1;
	}

	auto global_end = std::chrono::system_clock::now();
	chrono::duration<double> global_elapsed_seconds = global_end - global_start;
	cout << "Color elapsed time: " << global_elapsed_seconds.count() << "s\n";
	cout << endl;

	cout<<"kmer comparison calls: " << pretty_int(kmer_comp_call)<<endl;

	if (check)
		verif_colors(colors);

	cout << "Global statistics:" << endl;
	uint64_t nb_buckets, nb_skmers, nb_kmers, nb_cursed, memory, largest_bucket;
	colors.stats(nb_buckets, nb_skmers, nb_kmers, nb_cursed, memory, largest_bucket);
	cout << pretty_int(nb_buckets) << " bucket used (/" << pretty_int(pow(4, colors.params.m_small)) << " possible)" << endl;
	cout << "nb superkmers: " << pretty_int(nb_skmers) << endl;
	cout << "nb kmers: " << pretty_int(nb_kmers) << endl;
	cout << "kmer / second: " << pretty_int((float)nb_kmers / global_elapsed_seconds.count()) << endl;
	cout << "average kmer / superkmer: " << ((float)nb_kmers / (float)nb_skmers) << endl;
	cout << "average superkmer / bucket: " << ((float)nb_skmers / (float)nb_buckets) << endl;
	cout << "Largest bucket :	"<<pretty_int(largest_bucket) <<endl;
	cout << "Memory usage: " << (memory / 1024) << "Mo" << endl;
	cout << "bits / kmer: " << ((float)(memory * 1024 * 8) / (float)nb_kmers) << endl;
	cout << "nb cursed kmers: " << pretty_int(nb_cursed) << endl;
	cout<<"Nb kmer considered: " <<pretty_int(number_kmer_count)<<endl;
	cout<<"Low complexity kmer : "<<pretty_int(low_complexity_kmer)<<endl;
	cout<<"Comparison made : "<<pretty_int(comparaisons)<<endl;

	// --- Save Brisk index ---
	if (mode == 0 and outfile != "") {
		BriskWriter writer(outfile);
		writer.write(colors);
		writer.close();
	}

	return 0;
}


/** Counter function.
  * Read a complete fasta file line by line and store the counts into the Brisk datastructure.
  */
void add_color_from_fasta(Brisk<uint64_t> & colors, string & filename, const uint threads, uint64_t color) {
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
	cout << filename << " " << filename.length() << endl;
	zstr::ifstream in(filename);
	vector<string>  buffer;

	// uint idx = 0;
	//#pragma omp parallel num_threads(threads)
	{
		while (in.good() or not buffer.empty()) {
			string line;
			//#pragma omp critical
			{
				if(not buffer.empty()){
					line=buffer[buffer.size()-1];
					buffer.pop_back();
				}else{
					if (in.good()) {
						line = getLineFasta(&in);
					} else{
						line = "";
					}
				}
				if(line.size()>1000000000){
					buffer.push_back(line.substr(0,line.size()/2));
					buffer.push_back(line.substr(line.size()/2-colors.params.k+1));
					line="";
				}
			}

			if (line != "") {
				add_color_from_sequence(colors, line, color);
			}
		}
	}
}



void add_color_from_sequence(Brisk<uint64_t> & colors, string & sequence, uint64_t color) {
	// Line too short
	if (sequence.size() < colors.params.k){
		return;
	}

	vector<kmer_full> superkmer;

	SuperKmerEnumerator enumerator(sequence, colors.params.k, colors.params.m);
	kint minimizer = enumerator.next(superkmer);
	while (superkmer.size() > 0) {
		colors.protect_data(superkmer[0]);
		// Add the values
		if (check) {
			for (kmer_full & kmer : superkmer) {
				#pragma omp critical
				{
					if (verif.count(kmer.kmer_s) == 0){
						verif[kmer.kmer_s] = 0;
					}
					verif[kmer.kmer_s] |= color;
				}
			}
		}
		vector<bool> newly_inserted;
		vector<uint64_t*> vec(colors.insert_superkmer(superkmer,newly_inserted));
		for(uint i(0); i < vec.size();++i){
			uint64_t * data_pointer(vec[i]);
			if(newly_inserted[i]){
				(*data_pointer) = color;
			}else{
				(*data_pointer) |= color;
			}
		}
		colors.unprotect_data(superkmer[0]);
		// Next superkmer
		superkmer.clear();
		minimizer = enumerator.next(superkmer);
		if(minimizer==0){
			return;
		}
	}
	// cout<<"done"<<endl;
}





void verif_colors(Brisk<uint64_t> & colors) {
	cout << "--- Start counting verification ---" << endl;
	// cout<<verif.size() <<endl;
	// kint mini_mask = (1 << (2 * counter.m)) - 1;
	kmer_full kmer(0,0, colors.params.m, false);
	// Count 
	while (colors.next(kmer)) {
		if (verif.count(kmer.kmer_s) == 0) {
			cout << "pas dans verif weird"<<endl;
			verif[kmer.kmer_s] = 0;
		}else{
			
		}


		uint64_t * count = colors.get(kmer);
		if (count == NULL) {
			cout<<"NULL COUNT"<<endl;
			print_kmer(kmer.minimizer, colors.params.m); cout << endl;
			cout << (uint)kmer.minimizer << endl;
			print_kmer(kmer.kmer_s, colors.params.k); cout << endl;
			// cin.get();
		}else{
			verif[kmer.kmer_s] ^= *count;
		}

		kmer.kmer_s = 0;
	}

	// Summary print
	uint errors = 0;
	for (auto & it : verif) {
		if (it.second != 0) {
			errors += 1;
			if (it.second > 0) {
				cout << "missing "; print_kmer(it.first, colors.params.k); cout << " " << (uint)it.second << endl;
			} else {
				cout << "too many "; print_kmer(it.first, colors.params.k); cout << " " << (uint)(-it.second) << endl;
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

