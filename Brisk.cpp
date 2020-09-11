#include "strict_fstream.hpp"
#include "zstr.hpp"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <cstring>
#include <fstream>
//~ #include <gatbl/common.hpp>
//~ #include <gatbl/kmer.hpp>
#include <iostream>
#include <math.h>
#include <mutex>
#include <stdint.h>
#include <stdio.h>
#include <string>
#include <sys/resource.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <tmmintrin.h>
#include <unistd.h>
#include <unordered_map>
#include <vector>
#include "Kmers.cpp"
#include "SuperKmerLight.cpp"
#include "buckets.cpp"
#include "DenseMenuYo.hpp"
#include "pow2.hpp"
#include "robin_hood.h"



using namespace std;



DenseMenuYo menu(minimizer_size);
// SparseMenu menu(minimizer_size);
uint64_t line_count(0);
Pow2<kint> offsetUpdateAnchor(2 * k);
const Pow2<uint64_t> offsetUpdateAnchorMin(2 * super_minimizer_size);
// uint16_t abundance_mini[1<<(2*minimizer_size)];
// vector<Bucket> bucket_menus[1<<(2*subminimizer_size)];
uint64_t nb_core(20);



//START LOW LEVEL FUNCTIONS
uint64_t hash64shift(uint64_t key) {
	// uint64_t ab=abundance_mini[key];
	// uint64_t ab=0;
	// ab<<=32;
	key = (~key) + (key << 21); // key = (key << 21) - key - 1;
	key = key ^ (key >> 24);
	key = (key + (key << 3)) + (key << 8); // key * 265
	key = key ^ (key >> 14);
	key = (key + (key << 2)) + (key << 4); // key * 21
	key = key ^ (key >> 28);
	key = key + (key << 31);
	return (uint64_t)key;
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
	return result;
}



uint32_t revhash(uint32_t x) {
	x = ((x >> 16) ^ x) * 0x2c1b3c6d;
	x = ((x >> 16) ^ x) * 0x297a2d39;
	x = ((x >> 16) ^ x);
	return x;
}



uint32_t unrevhash(uint32_t x) {
	x = ((x >> 16) ^ x) * 0x0cf0b109; // PowerMod[0x297a2d39, -1, 2^32]
	x = ((x >> 16) ^ x) * 0x64ea2d65;
	x = ((x >> 16) ^ x);
	return x;
}



uint64_t revhash(uint64_t x) {
	x = ((x >> 32) ^ x) * 0xD6E8FEB86659FD93;
	x = ((x >> 32) ^ x) * 0xD6E8FEB86659FD93;
	x = ((x >> 32) ^ x);
	return x;
}



uint64_t unrevhash(uint64_t x) {
	x = ((x >> 32) ^ x) * 0xCFEE444D8B59A89B;
	x = ((x >> 32) ^ x) * 0xCFEE444D8B59A89B;
	x = ((x >> 32) ^ x);
	return x;
}



inline uint64_t nuc2int(char c) {
	return (c / 2) % 4;
}



inline uint64_t nuc2intrc(char c) {
	return ((c / 2) % 4) ^ 2;
}



inline void updateK(kint& min, char nuc) {
	min <<= 2;
	min += nuc2int(nuc);
	min %= offsetUpdateAnchor;
}



inline void updateM(uint64_t& min, char nuc) {
	min <<= 2;
	min += nuc2int(nuc);
	min %= offsetUpdateAnchorMin;
}



inline void updateRCK(kint& min, char nuc) {
	min >>= 2;
	min += ((kint)nuc2intrc(nuc) << (2 * k - 2));
}



inline void updateRCM(uint64_t& min, char nuc) {
	min >>= 2;
	min += (nuc2intrc(nuc) << (2 * super_minimizer_size - 2));
}


void clean(string& str){
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



void count_line(string& line) {
	//~ cout<<"COUNT LINE"<<endl;
	if (line.size() < k) {
		return;
	}
	clean(line);
	vector<kmer_full> kmers;
	//~ cout<<"count line"<<endl;
	// Init Sequences
	kint kmer_seq = (str2num(line.substr(0, k))), kmer_rc_seq(rcb(kmer_seq));
	uint64_t min_seq  = (uint64_t)(str2num(line.substr(k - super_minimizer_size, super_minimizer_size))), min_rcseq(rcbc(min_seq, super_minimizer_size)), min_canon(min(min_seq, min_rcseq));
	// Init MINIMIZER
	int8_t relative_min_position;
	int64_t minimizer = get_minimizer(kmer_seq, relative_min_position);
	bool multiple_min  = relative_min_position < 0;

	uint8_t position_minimizer_in_kmer;
	
	if (minimizer<0){
		position_minimizer_in_kmer = (uint8_t)(-relative_min_position - 1);
	}else{
		position_minimizer_in_kmer = relative_min_position;
	}

	uint64_t hash_mini = hash64shift(abs(minimizer));
	if (multiple_min) {
		cursed_kmers[min(kmer_rc_seq,kmer_seq)]++;
		//~ cout<<"CURSED"<<endl;
		//~ print_kmer(kmer_seq,31);cout<<endl;
	} else {
		if(minimizer<0){
			kmers.push_back({k-relative_min_position-super_minimizer_size+4, kmer_rc_seq});
		}else{
			if(kmer_seq!=0){//PUT COMPLEXITY THESHOLD
				kmers.push_back({relative_min_position+4, kmer_seq});
			}
		}
	}

	if (check) {
		real_count[getCanonical(line.substr(0, k))]++;
	}

	uint64_t line_size = line.size();
	for (uint64_t i = 0; i + k < line_size; ++i) {

		// Update KMER and MINIMIZER candidate with the new letter
		updateK(kmer_seq, line[i + k]);
		updateRCK(kmer_rc_seq, line[i + k]);
		updateM(min_seq, line[i + k]);
		updateRCM(min_rcseq, line[i + k]);
		min_canon = (min(min_seq, min_rcseq));

		//THE NEW mmer is a MINIMIZER
		uint64_t new_hash = (hash64shift(min_canon));
		//the previous MINIMIZER is outdated
		if (position_minimizer_in_kmer > k - super_minimizer_size) {
			cout << "Outdated mini" << endl;
			if(minimizer<0){
				reverse(kmers.begin(),kmers.end());
				minimizer*=-1;
			}
			menu.add_kmers(kmers,minimizer/256);
			// Search for the new MINIMIZER in the whole kmer
			minimizer    = get_minimizer(kmer_seq, relative_min_position);
			multiple_min = (relative_min_position < 0);
			if (multiple_min){
				position_minimizer_in_kmer = (uint8_t)(-relative_min_position );
			}else{
				position_minimizer_in_kmer = relative_min_position;
			}
			hash_mini = hash64shift(abs(minimizer));
		}
		else if (new_hash < hash_mini) {
			cout << "New mini" << endl;
			// Clear the previous kmer list
			if(minimizer<0){
				reverse(kmers.begin(),kmers.end());
				minimizer*=-1;
			}
			menu.add_kmers(kmers,minimizer/256);
			// Create the new minimizer
			minimizer                  = (min_canon);
			if(min_canon!=min_seq){minimizer*=-1;}
			hash_mini                  = new_hash;
			position_minimizer_in_kmer = relative_min_position = 0;
			multiple_min                                       = false;
		}
		// duplicated MINIMIZER
		else if (new_hash == hash_mini) {
			cout << "Duplicate mini" << endl;
			multiple_min = true;
			position_minimizer_in_kmer ++;
			relative_min_position = -((int8_t)position_minimizer_in_kmer) - 1;
			//~ cout<<relative_min_position<<endl;
		}
		else {
			//~ cout<<(int)position_minimizer_in_kmer<<" "<<(int)k - (int)super_minimizer_size-1<<endl;
			cout << "Nothing special" << endl;
			position_minimizer_in_kmer++;
			if (multiple_min){
				relative_min_position--;
			}else{
				relative_min_position++;
			}
		}

		// TODO: Multi-minimizer process
		if (multiple_min) {
			cursed_kmers[min(kmer_rc_seq,kmer_seq)]++;
			//~ cout<<"CURSED"<<endl;
			//~ print_kmer((uint64_t)kmer_rc_seq,31);cout<<endl;
		} else {
			// Normal add of the kmer into kmer list
			if(minimizer<0){
				int8_t val = ((int8_t)k) - ((int8_t)relative_min_position) - ((int8_t)super_minimizer_size) + 4;
				// kmers.push_back({k-relative_min_position-super_minimizer_size+4, kmer_rc_seq});
				kmers.push_back({val, kmer_rc_seq});
			}else{
				if(kmer_seq!=0){//PUT COMPLEXITY THESHOLD
					kmers.push_back({relative_min_position+4, kmer_seq});
				}
			}
		}
		if (check) {
			real_count[getCanonical(line.substr(i + 1, k))]++;
		}
	}
	if(minimizer<0){
		reverse(kmers.begin(),kmers.end());
		minimizer*=-1;
	}
	menu.add_kmers(kmers,minimizer/256);
	//~ menu.dump_counting();
}



inline bool exists_test (const string& name) {
  struct stat buffer;
  return (stat (name.c_str(), &buffer) == 0);
}



void read_fasta_file(const string& filename) {
	if(not exists_test(filename)){
		cout<<"Problem with file opening:	"<<filename<<endl;
		return;
	}
	zstr::ifstream in(filename);
	if (check) {
		nb_core = 1;
	}
	vector<string>  buffer;
	//#pragma omp parallel num_threads(nb_core)
	{
		string line;
		while (in.good() or  not buffer.empty()) {
			#pragma omp critical(input)
			{
				if(not buffer.empty()){
					line=buffer[buffer.size()-1];
					buffer.pop_back();
				}else{
					line = getLineFasta(&in);
					if(line.size()>100000000000){
						buffer.push_back(line.substr(0,line.size()/4));
						buffer.push_back(line.substr(line.size()/4-k+1,line.size()/4+k-1));
						buffer.push_back(line.substr(line.size()/2-k+1,line.size()/4+k-1));
						line=line.substr(3*line.size()/4-k+1);
					}
				}

			}
			count_line(line);
			line_count++;
		}
	}
}


	#include <bitset>

int main(int argc, char** argv) {
	if (argc < 2) {
		cout << "[fasta file]" << endl;
		cout << "0 DEFAULT, OUTPUT THE COUNT WITH NO CHECKING" << endl;
		cout << "1 PERFORMANCE MODE NO COUNT OUTPUT" << endl;
		cout << "2 DEBUG MODE COUNT AND CHECK WITH HASHTABLE" << endl;
		exit(0);
	}
	int mode(2);
	if (argc > 2) {
		mode = stoi(argv[2]);
	}

	if (mode > 1) {
		check = true;
		cout << "LETS CHECK THE RESULTS" << endl;
	}

	// SKCL skcl(0xFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFul, 0, 0);
	// cout << skcl.interleaved_value() << endl;
	// bitset<64> x(skcl.interleaved_value());
	// cout << x << endl;
	// return 0;


	cout << "\n\n\nI count " << argv[1] << endl;
	cout << "Minimizer size:	" << minimizer_size << endl;
	// cout << "Number of bucket:	" << bucket_number.value() << endl;
	auto start = std::chrono::system_clock::now();
	// read_fasta_file_ab(argv[1]);
	auto end = std::chrono::system_clock::now();
	chrono::duration<double> elapsed_seconds = end - start;
    // cout << "Minimize weight computed elapsed time: " << elapsed_seconds.count() << "s\n";

	start = std::chrono::system_clock::now();
	read_fasta_file(argv[1]);
	end = std::chrono::system_clock::now();
	elapsed_seconds = end - start;
	cout << "Kmer counted elapsed time: " << elapsed_seconds.count() << "s\n";
	cout << endl;
	if (mode % 2 == 0) {
		cout<<menu.dump_counting()<<" errors"<<endl;;
	}
	// cout<<"DUMP COUNTING DONE"<<endl;
	menu.dump_stats();
	exit(0);
}
