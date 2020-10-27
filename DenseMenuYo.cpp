// #include "kmers.hpp"
// #include "DenseMenuYo.hpp"

// DenseMenuYo::DenseMenuYo(uint8_t k, uint8_t m){
// 	// Create the mutexes
// 	mutex_order = min((uint8_t)6,m);
// 	mutex_number = 1<<(2*mutex_order); /* DO NOT MODIFY, to increase/decrease modify mutex_order*/
// 	// Init the mutexes
// 	MutexWall.resize(mutex_number);
// 	for (uint64_t i(0); i < mutex_number; ++i) {
// 		omp_init_lock(&MutexWall[i]);
// 	}

// 	minimizer_size=m;
// 	this->k = k;
// 	bucket_number=1<<(2*minimizer_size);
// 	matrix_column_number = bucket_number / mutex_number;
// 	// bucket_indexes = new uint32_t[mutex_number][bucket_number/mutex_number];
// 	bucket_indexes = new uint32_t[bucket_number];
// 	bucketMatrix = new vector<Bucket>[mutex_number];
// 	skm_total_size=call_ad=0;
// }


// #define get_mutex(mini) (mini%mutex_number)
// #define get_column(mini) (mini/mutex_number)
// #define matrix_position(row_idx, col_idx) (row_idx * matrix_column_number + col_idx)

// void DenseMenuYo::add_kmer(kmer_full & v, uint64_t minimizer){
// 	// // print_kmer(minimizer, minimizer_size);
// 	// // cout << " " << v.size() << endl;
// 	// #pragma omp critical
// 	// {
// 	// 	call_ad++;
// 	// 	skm_total_size+=v.size();
// 	// 	// size_sk[v.size()]++;
// 	// }
// 	// uint64_t mutex_idx = get_mutex(minimizer);
// 	// uint64_t column_idx = get_column(minimizer);
// 	// uint64_t matrix_idx = matrix_position(mutex_idx, column_idx);
// 	// omp_set_lock(&MutexWall[mutex_idx]);
// 	// uint32_t idx = bucket_indexes[matrix_idx];
// 	// if (idx == 0) {
// 	// 	bucketMatrix[mutex_idx].push_back(Bucket());
// 	// 	bucket_indexes[matrix_idx] = bucketMatrix[mutex_idx].size();
// 	// 	idx = bucket_indexes[matrix_idx];
// 	// }
// 	// bucketMatrix[mutex_idx][idx-1].add_kmers(v);
// 	// omp_unset_lock(&MutexWall[mutex_idx]);
// 	// v.clear();
// }

// void DenseMenuYo::add_kmers(vector<kmer_full>& v,uint64_t minimizer){
// 	// print_kmer(minimizer, minimizer_size);
// 	// cout << " " << v.size() << endl;
// 	#pragma omp critical
// 	{
// 		call_ad++;
// 		skm_total_size+=v.size();
// 		// size_sk[v.size()]++;
// 	}
// 	uint64_t mutex_idx = get_mutex(minimizer);
// 	uint64_t column_idx = get_column(minimizer);
// 	uint64_t matrix_idx = matrix_position(mutex_idx, column_idx);
// 	omp_set_lock(&MutexWall[mutex_idx]);
// 	uint32_t idx = bucket_indexes[matrix_idx];
// 	if (idx == 0) {
// 		bucketMatrix[mutex_idx].push_back(Bucket());
// 		bucket_indexes[matrix_idx] = bucketMatrix[mutex_idx].size();
// 		idx = bucket_indexes[matrix_idx];
// 	}
// 	bucketMatrix[mutex_idx][idx-1].add_kmers(v);
// 	omp_unset_lock(&MutexWall[mutex_idx]);
// 	v.clear();
// }



// uint64_t DenseMenuYo::dump_counting(bool check){
// 	uint64_t counting_errors = 0;

// 	string toprint;
// 	for(uint64_t mini(0);mini<bucket_number;++mini){
// 		uint64_t mutex_idx = get_mutex(mini);
// 		uint64_t column_idx = get_column(mini);
// 		uint64_t matrix_idx = matrix_position(mutex_idx, column_idx);
// 		uint32_t i = bucket_indexes[matrix_idx];
// 		if(i!=0){
// 			i--;
// 			if(bucketMatrix[mutex_idx][i].size()!=0){
// 				counting_errors += bucketMatrix[mutex_idx][i].print_kmers(
// 						toprint,
// 						kmer2str(mini,minimizer_size),
// 						real_counts,
// 						cursed_kmers,
// 						check
// 				);
// 			}
// 		}
// 	}
// 	if(check){
// 		for(int i(0);i<1024;++i){
// 			for (auto e:cursed_kmers[i]) {
// 				if(e.second!=0){
// 					string canonstr(getCanonical(kmer2str(e.first,k)));
// 					if(real_counts.count(canonstr)==1){
// 						if(real_counts[canonstr]==0){
// 							cout<<"This is a kmer in cursed AND in the main index "<<endl;
// 							cout<<canonstr<<" real "<<(int)real_counts[canonstr]<<" estimated:	"<<(int)e.second<<endl;
// 						}
// 					}
// 					if(real_counts[canonstr]!=e.second){
// 						cout<<"Error in cursed"<<endl;
// 						cout<<canonstr<<" real"<<(int)real_counts[canonstr]<<" estimated:	"<<(int)e.second<<endl;
// 						counting_errors++;
// 					}else{
// 						//~ cout<<"CURSED OK"<<endl;
// 					}
// 					real_counts[canonstr]=0;
// 				}
// 			}
// 		}
// 		for (auto e:real_counts) {
// 			if(e.second!=0){
// 				cout<<"I forgot	"<<e.first<<" "<<e.second<<endl;
// 				counting_errors++;
// 			}
// 		}
// 		if(counting_errors!=0){
// 			cout<< counting_errors<<"	errors"<<endl;
// 		}else{
// 			cout << "The results are OK" << endl;
// 		}
// 	}
// 	return counting_errors;
// }



// void DenseMenuYo::dump_stats(){
// 	uint64_t total_super_kmers(0);
// 	uint64_t total_kmers(0);
// 	uint64_t total_kmers_counted(0);
// 	uint64_t non_null_buckets(0);
// 	uint64_t null_buckets(0);
// 	uint64_t largest_bucket(0);
// 	for(uint64_t mini(0);mini<bucket_number;++mini){
// 		// cout<<"lini"<<mini<<endl;
// 		uint64_t mutex_idx = get_mutex(mini);
// 		uint64_t column_idx = get_column(mini);
// 		uint64_t matrix_idx = matrix_position(mutex_idx, column_idx);
// 		uint32_t i = bucket_indexes[matrix_idx];
// 		if(i!=0){
// 			i--;
// 			// cout<<"go"<<i<<endl;
// 			largest_bucket = max(largest_bucket,bucketMatrix[mutex_idx][i].size());
// 			// cout<<"lol1"<<endl;
// 			non_null_buckets++;
// 			total_super_kmers +=bucketMatrix[mutex_idx][i].size();
// 			// cout<<"lol2"<<endl;
// 			// cout<<bucketList[i].size()<<endl;
// 			total_kmers += bucketMatrix[mutex_idx][i].number_kmer();
// 			total_kmers_counted += bucketMatrix[mutex_idx][i].number_kmer_counted();
// 			// cout<<i<<"end"<<endl;
// 		}else{
// 			null_buckets++;
// 		}
// 	}
// 	cout << endl;
// 	cout << "Empty buckets:	" << intToString(null_buckets) << endl;
// 	cout << "Useful buckets:	" << intToString(non_null_buckets) << endl;
// 	cout << "#Superkmer:	" << intToString(total_super_kmers) << endl;
// 	// cout << "#Superkmer2:	" << intToString(nb_superkmer) << endl;
// 	cout << "#kmer:	" << intToString(total_kmers) << endl;
// 	cout << "#kmer (Counted):	" << intToString(total_kmers_counted) << endl;
// 	if(total_kmers!=0){
// 		cout << "super_kmer per useful buckets*1000:	" << intToString(total_super_kmers * 1000 / non_null_buckets) << endl;
// 		cout << "kmer per useful buckets*1000:	" << intToString(total_kmers * 1000 / non_null_buckets) << endl;
// 		cout << "kmer per super_kmer*1000:	" << intToString(total_kmers * 1000 / total_super_kmers) << endl;
// 		cout << "Largest_bucket:	" << intToString(largest_bucket) << endl;
// 		cout<<intToString(getMemorySelfMaxUsed())<<endl;
// 		cout<<intToString(getMemorySelfMaxUsed()*1024)<<" Bytes"<<endl;
// 		cout<<intToString(getMemorySelfMaxUsed()*1024*8/total_kmers)<<" Bits per kmer"<<endl;
// 		cout<<intToString(getMemorySelfMaxUsed()*1024/total_super_kmers)<<" Bytes per superkmer"<<endl;
// 		cout<<intToString(skm_total_size*1000/call_ad)<<" Real superkmer size"<<endl;
// 		uint64_t cursed_kmers_toal(0);
// 		for(int i(0);i<1024;++i){
// 			cursed_kmers_toal+=cursed_kmers[i].size();
// 		}
// 		cout<<"Number of cursed kmer	"<<intToString(cursed_kmers_toal)<<endl;
// 		// for (auto& it: size_sk) {
// 		// 	cout<<it.first<<" "<<intToString((uint64_t)it.second)<<endl;
// 		// }
// 	}
// 	cout<<"size of a super kmer:	"<<sizeof(SKCL)<<endl;
// 	// cout<<sizeof(uint256_t)<<endl;
// }


// uint64_t DenseMenuYo::getMemorySelfMaxUsed (){
// 	u_int64_t result = 0;
// 	struct rusage usage;
// 	if (getrusage(RUSAGE_SELF, &usage)==0)  {  result = usage.ru_maxrss;  }
// 	return result;
// }