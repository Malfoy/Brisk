#include "lib/kff/kff_io.hpp"
#include "Brisk.hpp"
#include "hashing.hpp"


#ifndef BRISKWRITER_H
#define BRISKWRITER_H

class BriskWriter {
private:
	Kff_file * current_file;
public:
	BriskWriter(std::string filename);
	template <class DATA>
	void write(Brisk<DATA> & index);
	void close();
};


BriskWriter::BriskWriter(std::string filename) {
	current_file = new Kff_file(filename, "w");
	// Set encoding   A  C  G  T
	current_file->write_encoding(0, 1, 3, 2);
	// Set metadata
	std::string meta = "File generated with Brisk v1. See https://github.com/Malfoy/Brisk";
	current_file->write_metadata(meta.length(), (uint8_t *)meta.c_str());
}

void little_to_big_endian(uint8_t * little, uint8_t * big, size_t bytes_to_convert) {
	for (uint8_t i=0 ; i<bytes_to_convert ; i++) {
		big[i] = little[bytes_to_convert - i - 1];
	}
}


/** Transform a pair of prefix suffix nucleotide uints into a big endian byte array.
  * If the number of nucleotides is not a multiple of 4, the bigest empty bits are set to 0.
 **/
void to_big_endian_compact(uint8_t * big, kint prefix, size_t pref_size, kint suffix, size_t suff_size) {
	uint64_t compact_size = pref_size + suff_size;
	uint64_t compact_bytes = (compact_size + 3) / 4;
	memset(big, 0, compact_bytes);

	uint64_t current_nucleotide = compact_bytes * 4 - 1;
	// Insert the suffix
	for (size_t rev_suff_idx=0 ; rev_suff_idx<suff_size ; rev_suff_idx++, current_nucleotide-- ) {
		size_t big_byte = current_nucleotide / 4;
		size_t big_nucl_pos = 3 - (current_nucleotide % 4);

		uint8_t nucl = suffix & 0b11;
		big[big_byte] |= nucl << (2 * big_nucl_pos);

		suffix >>= 2;
	}
	
	// Insert the prefix
	for (size_t rev_pref_idx=0 ; rev_pref_idx<pref_size ; rev_pref_idx++, current_nucleotide-- ) {
		size_t big_byte = current_nucleotide / 4;
		size_t big_nucl_pos = 3 - (current_nucleotide % 4);

		uint8_t nucl = prefix & 0b11;
		big[big_byte] |= nucl << (2 * big_nucl_pos);

		prefix >>= 2;
	}

	// cout << "compact ";
	// for (size_t i=0 ; i<compact_bytes ; i++)
	// 	cout << (uint64_t)big[i] << " ";
	// cout << endl;
}


static void rightshift8(uint8_t * bitarray, size_t length, size_t bitshift) {
	assert(bitshift < 8);

	if (length > 0) {
		for (uint64_t i=length-1 ; i>0 ; i--) {
			bitarray[i] = (bitarray[i-1] << (8-bitshift)) | (bitarray[i] >> bitshift);
		}
		bitarray[0] >>= bitshift;
	}
}


template <class DATA>
void BriskWriter::write(Brisk<DATA> & index) {
	// Set global variables
	Section_GV sgv(current_file);
	sgv.write_var("k", index.params.k);
	sgv.write_var("data_size", sizeof(DATA));
	sgv.write_var("max", 1);
	sgv.close();

	uint64_t nb_kmers = 0;

	// Save the cursed kmers into a raw block
	Section_Raw sr(current_file);
	DenseMenuYo<DATA> * menu = index.menu;
	uint8_t big_endian[32];
	uint8_t biggest_usefull_byte = index.params.k % 4 == 0 ? index.params.k / 4 : index.params.k / 4 + 1;
	for (auto it = menu->cursed_kmers.begin(); it != menu->cursed_kmers.end(); ++it) {
		kint kmer = it->first;
		// DATA & data = it.value();
		DATA & data = it->second;
		little_to_big_endian((uint8_t *)(&kmer), big_endian, biggest_usefull_byte);
		sr.write_compacted_sequence(big_endian, index.params.k, (uint8_t *)&data);
		nb_kmers += 1;
	}
	sr.close();

	// Prepare max value for super kmer size
	sgv = Section_GV(current_file);
	sgv.write_var("k", index.params.k);
	sgv.write_var("m", index.params.m);
	sgv.write_var("data_size", sizeof(DATA));
	sgv.write_var("max", 2 * (index.params.k - index.params.m));
	sgv.close();

	cout << "params: k=" << (uint64_t)index.params.k << " m=" << (uint64_t)index.params.m << " m_small=" << (uint64_t)index.params.m_small << " m_reduc=" << (uint64_t)index.params.m_reduc << endl;

	// Prepare structures for minimizer enumeration
	uint8_t bytes_mini = (index.params.m + 3) / 4;
	uint8_t * mini_seq = new uint8_t[bytes_mini];
	uint32_t max_mini = (1 << (2 * index.params.m_small)) - 1;
	// Enumerate all possible minimizers
	for (uint32_t small_mini=0 ; small_mini<=max_mini ; small_mini++) {
		uint32_t mutex_idx = small_mini % menu->mutex_number;
		uint32_t column_idx = small_mini >> (2 * menu->mutex_order);
		uint64_t matrix_idx = mutex_idx * menu->matrix_column_number + column_idx;
		uint32_t idx = menu->bucket_indexes[matrix_idx];

		// If the bucket for the minimizer exists
		if (idx != 0) {
			cout << "coordinates " << small_mini << " " << idx << endl;
			Section_Minimizer * sm = nullptr;
			uint64_t current_minimizer = 0xFFFFFFFFFFFFFFFFUL;

			uint8_t * big_endian_nucleotides = new uint8_t[index.params.allocated_bytes];
			Bucket<DATA> & b = menu->bucketMatrix[mutex_idx][idx-1];
			cout << "bucket " << b.skml.size() << endl;
			for (SKL & skmer : b.skml) {
				nb_kmers += skmer.size;
				// Get the right pointers
				uint8_t * nucleotides_ptr = b.nucleotides_reserved_memory + skmer.idx * index.params.allocated_bytes;
				DATA * data_ptr = b.data_reserved_memory + skmer.data_idx;

				// Get prefix and suffix
				kint prefix = skmer.get_prefix(nucleotides_ptr, index.params);
				kint suffix = skmer.get_suffix(nucleotides_ptr, index.params);

				// Minimizer prefix
				uint64_t mini_prefix_size = index.params.m_reduc / 2;
				uint64_t mask = (((kint)1) << (2 * mini_prefix_size)) - 1;
				uint64_t minimizer = prefix & mask;
				// Update prefix
				prefix >>= 2 * mini_prefix_size;
				// Small minimizer
				minimizer <<= index.params.m_small * 2;
				minimizer += small_mini;
				// Minimizer suffix
				uint64_t mini_suffix_size = (index.params.m_reduc + 1) / 2;
				mask = (((kint)1) << (2 * mini_suffix_size)) - 1;
				minimizer <<= mini_suffix_size * 2;
				cout << (skmer.suffix_size() - mini_suffix_size) << endl;
				minimizer += (suffix >> (2 * (skmer.suffix_size() - mini_suffix_size))) & mask;
				// Update suffix
				mask = (((kint)1) << (2 * (skmer.suffix_size() - mini_suffix_size))) - 1;
				suffix &= mask;
				// Unhash minimizer
				minimizer = bfc_hash_64_inv(minimizer, index.params.mask_large_minimizer);

				cout << "skmer " << (uint64_t)skmer.size << endl;
				cout << "sizes: " << "prefix=" << (uint64_t)skmer.prefix_size(index.params);
				cout << " mini_prefix_size=" << (uint64_t)mini_prefix_size;
				cout << " m=" << (uint64_t)index.params.m;
				cout << " suffix=" << (uint64_t)skmer.suffix_size();
				cout << " mini_suffix_size=" << (uint64_t)mini_suffix_size << endl;
				cout << kmer2str(prefix, skmer.prefix_size(index.params) - mini_prefix_size);
				cout << " " << kmer2str(minimizer, index.params.m) << " ";
				cout << kmer2str(suffix, skmer.suffix_size() - mini_suffix_size) << endl;

				// Create a new minimizer section while changing
				if (current_minimizer != minimizer) {
					// Init new section object
					if (sm != nullptr) {
						sm->close();
						delete sm;
					}
					sm = new Section_Minimizer(current_file);
					memset(mini_seq, 0, bytes_mini);

					// Create the minimizer in bigendian order
					size_t pos = bytes_mini * 4 - 1;
					for (size_t nucl_idx=0 ; nucl_idx<index.params.m ; nucl_idx++, pos--) {
						size_t byte_idx = pos / 4;
						size_t byte_position = 3 - (pos % 4);
						uint8_t nucl = (minimizer >> (2 * nucl_idx)) & 0b11;
						mini_seq[byte_idx] |= nucl << (2 * byte_position);
					}

					// Write the minimizer inside of the section
					sm->write_minimizer(mini_seq);
				}

				// Transform the superkmer into a compacted big endian byte array
				to_big_endian_compact(big_endian_nucleotides,
									  prefix, skmer.prefix_size(index.params) - mini_prefix_size,
									  suffix, skmer.suffix_size() - mini_suffix_size);
				size_t real_seq_size = index.params.k + skmer.size - 1 - index.params.m;

				// Save the whole skmer at once
				sm->write_compacted_sequence_without_mini(
						big_endian_nucleotides,
						real_seq_size,
						skmer.prefix_size(index.params),
						// /!\ WARNING: This piece of code car trigger problems of little/big endian
						// We have to report it in the documentation
						(uint8_t *)data_ptr
				);
			}

			sm->close();
			delete sm;
			delete[] big_endian_nucleotides;
		}
	}

	delete[] mini_seq;
}





void BriskWriter::close() {
	current_file->close();
	delete current_file;
}

#endif
