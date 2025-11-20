#ifndef SEARCH_H
#define SEARCH_H

#include <vector>
#include <string>
#include "fasta.h"
#include "index.h"
#include "scoring.h"

// High Scoring Pair (HSP) structure
struct HSP {
    int sid;           // Sequence index in database
    int db_start;      // Start position in database
    int db_end;        // End position in database
    int q_start;       // Start position in query
    int q_end;         // End position in query
    int score;         // Alignment score
    double identity;   // Percent identity
};

// Find all HSPs for a query sequence
std::vector<HSP> findHSPs(
    const std::string& query,
    const std::vector<Sequence>& database,
    const KmerIndex& index,
    int k
);

// Merge overlapping HSPs for the same sequence
// Keeps the best scoring HSP when overlaps occur
std::vector<HSP> mergeHSPs(const std::vector<HSP>& hsps);

// Get alignment string representation
std::string getAlignment(
    const std::string& db_seq,
    const std::string& query,
    int db_start,
    int db_end,
    int q_start,
    int q_end
);

#endif // SEARCH_H

