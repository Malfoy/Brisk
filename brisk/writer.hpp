#include "lib/kff/kff_io.hpp"
#include "Brisk.hpp"


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

template <class DATA>
void BriskWriter::write(Brisk<DATA> & index) {
	// Set global variables
	Section_GV sgv(current_file);
	sgv.write_var("k", index.params.k);
	sgv.write_var("m", index.params.m_small);
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
	sgv.write_var("max", 2 * (index.params.k - index.params.m_small));
	sgv.close();

	// Prepare structures for minimizer enumeration
	uint8_t bytes_mini = index.params.m_small % 4 == 0 ? index.params.m_small / 4 : index.params.m_small / 4 + 1;
	uint8_t * mini_seq = new uint8_t[bytes_mini];
	uint32_t max_mini = (1 << (2 * index.params.m_small)) - 1;
	// Enumerate all possible minimizers
	for (uint32_t minimizer=0 ; minimizer<=max_mini ; minimizer++) {
		uint32_t mutex_idx = minimizer % menu->mutex_number;
		uint32_t column_idx = minimizer >> (2 * menu->mutex_order);
		uint64_t matrix_idx = mutex_idx * menu->matrix_column_number + column_idx;
		uint32_t idx = menu->bucket_indexes[matrix_idx];

		// If the bucket for the minimizer exists
		if (idx != 0) {
			Section_Minimizer sm(current_file);
			little_to_big_endian((uint8_t *)(&minimizer), mini_seq, bytes_mini);
			sm.write_minimizer(mini_seq);

			uint8_t * big_endian_nucleotides = new uint8_t[index.params.allocated_bytes];
			Bucket<DATA> & b = menu->bucketMatrix[mutex_idx][idx-1];
			for (SKL & skmer : b.skml) {
				nb_kmers += skmer.size;
				// Get the right pointers
				uint8_t * nucleotides_ptr = b.nucleotides_reserved_memory + skmer.idx * index.params.allocated_bytes;
				// TODO: Maybe a bug here
				DATA * data_ptr = b.data_reserved_memory + skmer.data_idx;

				// Little endian to big endian
				size_t real_seq_size = index.params.k + skmer.size - 1 - index.params.m_small;
				size_t real_seq_bytes = real_seq_size % 4 == 0 ? real_seq_size / 4 : real_seq_size / 4 + 1;
				// cout << real_seq_size << " " << real_seq_bytes << endl;

				// Copy the superkmer
				for (uint i=0 ; i<real_seq_bytes ; i++) {
					big_endian_nucleotides[i] = nucleotides_ptr[index.params.allocated_bytes-i-1];
				}
				// Align the skmer on the right
				if (real_seq_size % 4 != 0) {
					rightshift8(big_endian_nucleotides, real_seq_bytes, ((4 - (real_seq_size % 4)) % 4) * 2);
				}

				// Save the whole skmer at once
				sm.write_compacted_sequence_without_mini(
						big_endian_nucleotides,
						real_seq_size,
						skmer.prefix_size(index.params),
						data_ptr // TODO: Little/big endian
				);
			}

			sm.close();
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
