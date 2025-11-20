#ifndef SCORING_H
#define SCORING_H

#include <string>

// Ungapped extension result
struct ExtensionResult {
    int db_start;      // Start position in database sequence
    int db_end;        // End position in database sequence
    int q_start;       // Start position in query sequence
    int q_end;         // End position in query sequence
    int score;         // Alignment score
    double identity;   // Percent identity (0-100)
};

// Perform ungapped extension from a seed position
// Match: +2, Mismatch: -1
// Extends both left and right from seed
ExtensionResult extendUngapped(
    const std::string& db_seq,
    const std::string& query,
    int db_seed_pos,
    int q_seed_pos
);

// Calculate percent identity for an alignment
double calculateIdentity(
    const std::string& db_seq,
    const std::string& query,
    int db_start,
    int db_end,
    int q_start,
    int q_end
);

#endif // SCORING_H

