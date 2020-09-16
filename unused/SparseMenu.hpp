#include <iostream>
#include <cstdint>
#include "Kmers.hpp"
#include "SuperKmerCount.hpp"
#include "buckets.hpp"
#include "sparse_map.h"
#include <omp.h>





#ifndef SparseMENU_H
#define SparseMENU_H



const  uint64_t number_mutex(1024);



class SparseMenu{
public:
	robin_hood::unordered_flat_map<uint32_t,Bucket> min2b[number_mutex];
	// tsl::sparse_pg_map<uint64_t,Bucket> min2b[1024];
	uint64_t minimizer_size;
	omp_lock_t MutexWall[number_mutex];

	int dump_counting(){
		string toprint;
		for(uint i(0);i<number_mutex;++i){
			for (auto& B: min2b[i]) {
				B.second.print_kmers(toprint,kmer2str((uint64_t)B.first,minimizer_size));
			}
		}

		if(check){
			for (auto e:real_count) {
				if(e.second!=0){
					cout<<"I forgot	"<<e.first<<" "<<e.second<<endl;
					counting_errors++;
				}
			}
			if(counting_errors!=0){
				cout<< counting_errors<<"	errors"<<endl;
			}else{
				cout << "The results are OK" << endl;
			}
		}
		return counting_errors;
	}



	void dump_stats(){
		uint64_t total_super_kmers(0);
		uint64_t total_kmers(0);
		uint64_t non_null_buckets(0);
		uint64_t null_buckets(0);
		uint64_t largest_bucket(0);
		for(uint i(0);i<number_mutex;++i){
			for (auto& B: min2b[i]) {
				largest_bucket = max(largest_bucket,B.second.size());
				non_null_buckets++;
				total_super_kmers +=B.second.size();
				total_kmers +=B.second.number_kmer();
			}
		}
		cout << endl;
		cout << "Empty buckets:	" << intToString(null_buckets) << endl;
		cout << "Useful buckets:	" << intToString(non_null_buckets) << endl;
		cout << "#Superkmer:	" << intToString(total_super_kmers) << endl;
		// cout << "#Superkmer2:	" << intToString(nb_superkmer) << endl;
		cout << "#kmer:	" << intToString(total_kmers) << endl;
		cout << "super_kmer per useful buckets*1000:	" << intToString(total_super_kmers * 1000 / non_null_buckets) << endl;
		cout << "kmer per useful buckets*1000:	" << intToString(total_kmers * 1000 / non_null_buckets) << endl;
		cout << "kmer per super_kmer*1000:	" << intToString(total_kmers * 1000 / total_super_kmers) << endl;
		cout << "Largest_bucket:	" << intToString(largest_bucket) << endl;
		cout<<intToString(getMemorySelfMaxUsed()*1024)<<" Bytes"<<endl;
		cout<<intToString(getMemorySelfMaxUsed()*1024/total_kmers)<<" Bytes per kmer"<<endl;
		cout<<intToString(getMemorySelfMaxUsed()*1024/total_super_kmers)<<" Bytes per superkmer"<<endl;
		// cout<<sizeof(vector<Bucket>)<<endl;
	}



	void add_kmers(vector<kmer_full>& v,uint64_t minimizer){
		omp_set_lock(&MutexWall[minimizer % number_mutex]);
		//TODO WE CAN AVOID TO ENCODE THE  first log(bucket_number) bits
		min2b[minimizer % number_mutex][minimizer].add_kmers(v);
		// min2b[minimizer % 1024][minimizer]={};
		omp_unset_lock(&MutexWall[minimizer % number_mutex]);
		v.clear();
		// string toprint;
		// bucketList[minimizer].print_kmers(toprint,kmer2str(minimizer,minimizer_size));
		// if(dump_counting()!=0){
		// 	cout<<"EPIC FAIL"<<endl;
		// 	exit(0);
		// }
		// cout<<"ADDKMEREND"<<endl;
	}



	SparseMenu(uint64_t minisize){
		for (uint64_t i(0); i < number_mutex; ++i) {
			omp_init_lock(&MutexWall[i]);
		}
		minimizer_size=minisize;
	}
	uint64_t getMemorySelfMaxUsed (){
		u_int64_t result = 0;
		struct rusage usage;
		if (getrusage(RUSAGE_SELF, &usage)==0)  {  result = usage.ru_maxrss;  }
		return result;
	}
};

#endif
