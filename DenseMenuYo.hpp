#include <iostream>
#include <cstdint>
#include <cstring>
#include <map>
#include "Kmers.hpp"
#include "SuperKmerCount.hpp"
#include "buckets.hpp"
#include <omp.h>
#include <vector>





#ifndef DENSEMENUYO_H
#define DENSEMENUYO_H


#define mutex_order (6)
#define mutex_number (1<<(2*mutex_order)) /* DO NOT MODIFY, to increase/decrease modify mutex_order*/



class DenseMenuYo{
public:
	uint32_t * indexes;
	// std::vector<Bucket> bucketList;
	std::vector<Bucket> * bucketMatrix;
	uint64_t minimizer_size;
	uint64_t bucket_number;
	uint64_t matrix_column_number;
	omp_lock_t MutexWall[mutex_number];
	u_int64_t call_ad;
	u_int64_t skm_total_size;
	map<uint,uint> size_sk;


	DenseMenuYo(uint64_t minisize){
		for (uint64_t i(0); i < mutex_number; ++i) {
			omp_init_lock(&MutexWall[i]);
		}
		minimizer_size=minisize;
		bucket_number=1<<(2*minimizer_size);
		matrix_column_number = bucket_number / mutex_number;
		// indexes = new uint32_t[mutex_number][bucket_number/mutex_number];
		indexes = new uint32_t[bucket_number];
		bucketMatrix = new std::vector<Bucket>[mutex_number];
		skm_total_size=call_ad=0;
	}

	#define get_mutex(mini) (mini%mutex_number)
	#define get_column(mini) (mini/mutex_number)
	#define matrix_position(row_idx, col_idx) (row_idx * matrix_column_number + col_idx)

	void add_kmers(vector<kmer_full>& v,uint64_t minimizer){
		#pragma omp critical
		{
			call_ad++;
			skm_total_size+=v.size();
			size_sk[v.size()]++;
		}
		// cout << "add" << endl;
		uint64_t mutex_idx = get_mutex(minimizer);
		uint64_t column_idx = get_column(minimizer);
		uint64_t matrix_idx = matrix_position(mutex_idx, column_idx);
		omp_set_lock(&MutexWall[mutex_idx]);
			uint32_t idx = indexes[matrix_idx];
			// cout << idx << " " << (bucket_number/mutex_number) << " " << bucketMatrix[mutex_idx].size() << endl;
			if (idx == 0) {
				// cout << "new bucket" << endl;
				bucketMatrix[mutex_idx].push_back(Bucket());
				indexes[matrix_idx] = bucketMatrix[mutex_idx].size();
				idx = indexes[matrix_idx];
				// cout << "/new bucket" << endl;
			}
			// else cout << "fill bucket" << endl;
		// cout << "True add" << endl;
		// cout << idx << " " << (bucket_number/mutex_number) << " " << bucketMatrix[mutex_idx].size() << endl;
		bucketMatrix[mutex_idx][idx-1].add_kmers(v);
		omp_unset_lock(&MutexWall[mutex_idx]);

		v.clear();
		// cout << "/add" << endl;
	}


	int dump_counting(){
		string toprint;
		for(uint64_t mini(0);mini<bucket_number;++mini){
			uint64_t mutex_idx = get_mutex(mini);
			uint64_t column_idx = get_column(mini);
			uint64_t matrix_idx = matrix_position(mutex_idx, column_idx);
			uint32_t i = indexes[matrix_idx];
			if(i!=0){
				i--;
				if(bucketMatrix[mutex_idx][i].size()!=0){
					bucketMatrix[mutex_idx][i].print_kmers(toprint,kmer2str(mini,minimizer_size));
				}
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
		for(uint64_t mini(0);mini<bucket_number;++mini){
			// cout<<"lini"<<mini<<endl;
			uint64_t mutex_idx = get_mutex(mini);
			uint64_t column_idx = get_column(mini);
			uint64_t matrix_idx = matrix_position(mutex_idx, column_idx);
			uint32_t i = indexes[matrix_idx];
			if(i!=0){
				i--;
				// cout<<"go"<<i<<endl;
				largest_bucket = max(largest_bucket,bucketMatrix[mutex_idx][i].size());
				// cout<<"lol1"<<endl;
				non_null_buckets++;
				total_super_kmers +=bucketMatrix[mutex_idx][i].size();
				// cout<<"lol2"<<endl;
				// cout<<bucketList[i].size()<<endl;
				total_kmers += bucketMatrix[mutex_idx][i].number_kmer();
				// cout<<i<<"end"<<endl;
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
			cout<<intToString(skm_total_size*1000/call_ad)<<" Real superkmer size"<<endl;
			// for (auto& it: size_sk) {
			// 	cout<<it.first<<" "<<intToString((uint64_t)it.second)<<endl;
			// }
		}
		cout<<sizeof(SKCL)<<endl;
		// cout<<sizeof(uint256_t)<<endl;
	}


	uint64_t getMemorySelfMaxUsed (){
	    u_int64_t result = 0;
	    struct rusage usage;
	    if (getrusage(RUSAGE_SELF, &usage)==0)  {  result = usage.ru_maxrss;  }
	    return result;
	}




};

#endif
