#ifndef INDEX_H
#define INDEX_H

#include <unordered_map>
#include <vector>
#include <string>
#include "fasta.h"

// Hash table: k-mer key -> list of (sequence_index, position)
// Using uint32_t for k-mer encoding (supports k up to 16)
using KmerIndex = std::unordered_map<uint32_t, std::vector<std::pair<int, int>>>;

// Build k-mer hash index from database sequences
// Uses 2-bit encoding: A=0, C=1, G=2, T=3
// Uses rolling hash for efficient k-mer extraction
KmerIndex buildIndex(const std::vector<Sequence>& database, int k);

// Encode a k-mer string to integer using 2-bit encoding
// Each nucleotide takes 2 bits: A=00, C=01, G=10, T=11
uint32_t encodeKmer(const std::string& kmer);

// Extract k-mer at position i from sequence using rolling hash
// Returns encoded k-mer value
uint32_t getKmerAt(const std::string& seq, int pos, int k);

#endif // INDEX_H

