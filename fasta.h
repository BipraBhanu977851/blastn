#ifndef FASTA_H
#define FASTA_H

#include <string>
#include <vector>

// Structure to hold a database sequence with its metadata
struct Sequence {
    std::string id;           // Sequence ID (before |)
    std::string species;      // Species name (after |)
    std::string seq;          // DNA sequence
    int index;                // Index in database vector
};

// Structure to hold a query sequence with its name
struct Query {
    std::string name;         // Query name (from header)
    std::string seq;          // DNA sequence
};

// Parse database FASTA file with multiple sequences
// Format: >id|species\nsequence
std::vector<Sequence> parseDatabase(const std::string& filename);

// Parse query FASTA file (single sequence)
// Returns the DNA sequence string
std::string parseQuery(const std::string& filename);

// Parse query FASTA file with multiple sequences
// Returns vector of Query structures
std::vector<Query> parseQueries(const std::string& filename);

#endif // FASTA_H

