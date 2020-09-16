#include <iostream>
#include <cstdint>
#include "Kmers.hpp"
#include "SuperKmerCount.hpp"
#include "buckets.hpp"
#include <omp.h>





#ifndef DENSEMENU_H
#define DENSEMENU_H


class DenseMenu{
public:
	Bucket* bucketList;
	uint64_t minimizer_size;
	uint64_t bucket_number;
	omp_lock_t MutexWall[8192];


	int dump_counting(){
		string toprint;
		for(uint64_t i(0);i<bucket_number;++i){
			if(bucketList[i].size()!=0){
				bucketList[i].print_kmers(toprint,kmer2str(i,minimizer_size));
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
		for(uint64_t i(0);i<bucket_number;++i){
			if(bucketList[i].size()!=0){
				largest_bucket = max(largest_bucket,bucketList[i].size());
				non_null_buckets++;
				total_super_kmers +=bucketList[i].size();
				total_kmers += bucketList[i].number_kmer();
			}else{
				null_buckets++;
			}
		}
		cout << endl;
		cout << "Empty buckets:	" << intToString(null_buckets) << endl;
		cout << "Useful buckets:	" << intToString(non_null_buckets) << endl;
		cout << "#Superkmer:	" << intToString(total_super_kmers) << endl;
		// cout << "#Superkmer2:	" << intToString(nb_superkmer) << endl;
		cout << "#kmer:	" << intToString(total_kmers) << endl;
		if(total_kmers!=0){
			cout << "super_kmer per useful buckets*1000:	" << intToString(total_super_kmers * 1000 / non_null_buckets) << endl;
			cout << "kmer per useful buckets*1000:	" << intToString(total_kmers * 1000 / non_null_buckets) << endl;
			cout << "kmer per super_kmer*1000:	" << intToString(total_kmers * 1000 / total_super_kmers) << endl;
			cout << "Largest_bucket:	" << intToString(largest_bucket) << endl;
			cout<<intToString(getMemorySelfMaxUsed())<<endl;
			cout<<intToString(getMemorySelfMaxUsed()*1024)<<" Bytes"<<endl;
			cout<<intToString(getMemorySelfMaxUsed()*1024*8/total_kmers)<<" Bits per kmer"<<endl;
			cout<<intToString(getMemorySelfMaxUsed()*1024/total_super_kmers)<<" Bytes per superkmer"<<endl;
		}
		cout<<sizeof(SKC)<<endl;
	}



	void add_kmers(vector<kmer_full>& v,uint64_t minimizer){
		omp_set_lock(&MutexWall[minimizer % 8192]);
		bucketList[minimizer].add_kmers(v);
		// bucketList[minimizer]={};
		omp_unset_lock(&MutexWall[minimizer % 8192]);
		v.clear();
	}



	DenseMenu(uint64_t minisize){
		for (uint64_t i(0); i < 8192; ++i) {
			omp_init_lock(&MutexWall[i]);
		}
		minimizer_size=minisize;
		bucket_number=1<<(2*minimizer_size);
		bucketList=new Bucket[bucket_number];
	}



	uint64_t getMemorySelfMaxUsed (){
	    u_int64_t result = 0;
	    struct rusage usage;
	    if (getrusage(RUSAGE_SELF, &usage)==0)  {  result = usage.ru_maxrss;  }
	    return result;
	}
};

#endif
