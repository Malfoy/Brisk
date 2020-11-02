# Brisk - Brisk Reduced Index for Sequence of K-mers

Brisk is library for dynamic kmer indexing.
You can associate a data to each kmer, access it and modify it easily.
The kmer index can be used and update at the same time.

The library take adventage of overlapping kmers to store them in a compact way.

**Limitations**: The maximum k size is 64 due to underlying uint128 operations.

## Brisk datastructure overview overview

TODO

## Brisk API

### Create a Brisk index

A Brisk index is creating using 3 parameters:
- k: kmer length (up to 64)
- m: minimizer size (must be odd)
- b: bucket order of magniture (4^b minimizer / bucket)

```cpp
  #include "Brisk.hpp"

  uint8_t k = 63;
  uint8_t m = 11;
  uint8_t b = 4;

  Parameters params(k, m, b);

  // Create a Brisk index where the data associated to each kmer are uint8_t values.
  Brisk<uint8_t> index(params);
```

### Get/Insert kmers

There are two main operations that can be performed on a Brisk datastructure:
- Insert a kmer: If the kmer is absent from the index, allocate a DATA space corresponding to the kmer and return a pointer to it.  If the kmer is already present, similar to get function.
- Get a kmer data: Return the pointer to the coresponding kmer DATA. NULL if the kmer is absent from the index.

```cpp
  Brisk<uint8_t> index(params);
  kmer_full kmer = str2kmer("CTTAAAGAGATTTGCGGTCAACCGTTTTTTGAAAAAATTTTATAAAAATATTTATCATATTGT", params);

  // Get an absent kmer
  uint8_t * data_pointer = index.get(kmer);
  if (data_pointer == NULL)
    std::cout << "Absent kmer" << std::endl;
  
  // Insert the kmer and modify the value associated
  data_pointer = index.insert(kmer);
  *data_pointer = 21;

  // Get the data related to the kmer and print the value
  data_pointer = index.get(kmer);
  std::cout << (uint)*data_pointer << std::endl;
  // The value can also be updated after a get
  *data_pointer *= 2;
```

### Make it parallel (openMP)

## Example - Brisk counter index


