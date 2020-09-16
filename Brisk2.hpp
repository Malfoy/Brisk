#include <stdint.h>
#include <string>

class Brisk {
public:
	Brisk(uint8_t k, uint8_t m);

	void insert_sequence(std::string seq);
	void insert_fasta_file(std::string filename);
};
