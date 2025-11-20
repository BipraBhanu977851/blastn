#include "index.h"
#include <iostream>

// Encode a k-mer string to integer using 2-bit encoding
// Each nucleotide takes 2 bits: A=00, C=01, G=10, T=11
// This allows k up to 16 (32 bits / 2 bits per nucleotide)
uint32_t encodeKmer(const std::string& kmer) {
    uint32_t encoded = 0;
    
    for (char c : kmer) {
        encoded <<= 2;  // Shift left by 2 bits
        
        switch (c) {
            case 'A': case 'a':
                encoded |= 0;  // 00
                break;
            case 'C': case 'c':
                encoded |= 1;  // 01
                break;
            case 'G': case 'g':
                encoded |= 2;  // 10
                break;
            case 'T': case 't':
                encoded |= 3;  // 11
                break;
            default:
                // Invalid character, return 0 or handle error
                return 0;
        }
    }
    
    return encoded;
}

// Extract k-mer at position i from sequence using rolling hash
// This function encodes the k-mer starting at position pos
uint32_t getKmerAt(const std::string& seq, int pos, int k) {
    if (pos + k > static_cast<int>(seq.length())) {
        return 0;  // Invalid position
    }
    
    std::string kmer = seq.substr(pos, k);
    return encodeKmer(kmer);
}

// Build k-mer hash index from database sequences
// Uses 2-bit encoding: A=0, C=1, G=2, T=3
// Uses rolling hash for efficient k-mer extraction
KmerIndex buildIndex(const std::vector<Sequence>& database, int k) {
    KmerIndex index;
    
    // For each sequence in database
    for (const auto& seq : database) {
        const std::string& sequence = seq.seq;
        
        // Extract all k-mers using rolling window
        // For each position from 0 to seq.length() - k
        for (int i = 0; i <= static_cast<int>(sequence.length()) - k; ++i) {
            // Get k-mer at position i
            uint32_t kmer_key = getKmerAt(sequence, i, k);
            
            // Skip invalid k-mers (containing N or other invalid chars)
            if (kmer_key == 0 && k > 0) {
                // Check if it's actually invalid (not just a valid k-mer that encodes to 0)
                bool valid = true;
                for (int j = 0; j < k; ++j) {
                    char c = sequence[i + j];
                    if (c != 'A' && c != 'C' && c != 'G' && c != 'T' &&
                        c != 'a' && c != 'c' && c != 'g' && c != 't') {
                        valid = false;
                        break;
                    }
                }
                if (!valid) continue;
            }
            
            // Store (sequence_index, position) in hash table
            index[kmer_key].push_back({seq.index, i});
        }
    }
    
    return index;
}

