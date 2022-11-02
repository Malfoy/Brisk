#include <cstdint>


#ifndef HASHING_HPP
#define HASHING_HPP


uint64_t bfc_hash_64(uint64_t key, uint64_t mask);
uint64_t bfc_hash_64_inv(uint64_t key, uint64_t mask);


#endif