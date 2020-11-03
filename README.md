# Brisk - Brisk Reduced Index for Sequence of K-mers

Brisk is library for dynamic kmer indexing.
You can associate a data to each kmer, access it and modify it easily.
The kmer index can be used and update at the same time.

The library take advantage of overlapping kmers to store them in a compact way.

**Limitations**: The maximum k size is 64 due to underlying uint128 operations.

## Brisk datastructure overview

TODO

## Brisk API

### Create a Brisk index

A Brisk index is creating using 3 parameters:
- k: kmer length (up to 64)
- m: minimizer size (must be odd)
- b: bucket order of magnitude (4^b minimizer / bucket)

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
- Get a kmer data: Return the pointer to the corresponding kmer DATA. NULL if the kmer is absent from the index.

```cpp
  Brisk<uint8_t> index(params);
  // Warning: This way to create a kmer is not efficient. Always prefer to enumerate kmers from a sequence as in the next part.
  kmer_full kmer = str2kmer("CTTAAAGAGATTTGCGGTCAACCGTTTTTTGAAAAAATTTTATAAAAATATTTATCATATTGT", params.m);

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

### Enumerate kmers in sequence(s)

Kmers are stored as superkmers inside of the datastructure.
So, Brisk is very efficient when kmers are added by overlapping successors.
To help users to enumerate kmers in the right order we provide an enumerator object called SuperKmerEnumerator.
This object is able to translate a string into a succession of kmer vectors.

```cpp
  string line = readline_fasta(my_file_pointer); // Use here your favorite reader

  // Create the enumerator
  SuperKmerEnumerator enumerator(line, params.k, params.m);

  vector<kmer_full> superkmer;
  // The minimizer is a unsigned interger of kint type. for now kint is defined as a uint128_t
  auto minimizer = enumerator.next(superkmer);
  while (superkmer.size() > 0) {
    for (kmer_full & kmer : superkmer) {
      data_pointer = counter.insert(kmer);
      *data_pointer = 0;
    }

    // Next superkmer computation
    superkmer.clear();
    minimizer = enumerator.next(superkmer);
  }
```

If you don't want to use our enumerator, you can directly fill kmer_full objects.
But be careful when you insert kmers into Brisk: insert overlapping kmers will optimize space where random kmers will highly decrease space and time performances.

### Make it parallel (openMP)

The Brisk core code is openMP thread safe.
So get and insert can be done in parallel using pragma omp parallel (see counter.cpp for a full example).

But data pointers can be invalidate if insertions are performed by other threads and if chunks of memory are reallocated.
So, to prevent reallocation during data modifications you can protect and unprotect data as follow

```cpp
  // Run the code on 8 threads using openmp
  #pragma omp parallel num_threads(8)
  {
    while (not my_file_pointer.eof()) {
      string line = "";

      // Protect from file reading conflicts
      #pragma omp critical
      {
        line = readline_fasta(my_file_pointer); // Use here your favorite reader
      }

      // Create the enumerator
      SuperKmerEnumerator enumerator(line, params.k, params.m);
      vector<kmer_full> superkmer;
      // The minimizer is a unsigned interger of kint type. for now kint is defined as a uint128_t
      auto minimizer = enumerator.next(superkmer);
      while (superkmer.size() > 0) {
        for (kmer_full & kmer : superkmer) {
          // Forbid reallocations inside of the kmer bucket
          index.protect_data(kmer);
          
          // Get the data pointer and safely modify the value
          uint8_t * data_pointer = index.insert(kmer);
          *data_pointer = 84;

          // Release constraints on the kmer bucket
          index.unprotect_data(kmer);
        }

        // Next superkmer computation
        superkmer.clear();
        minimizer = enumerator.next(superkmer);
      }
    }
  }
```

**WARNING**: Different kmers with different minimizers can share the same bucket.
So, their protection can be mutually exclusive.
The protection of multiple kmers simultaneously in the same thread can provoke interlocking.

## Example - Brisk counter index


