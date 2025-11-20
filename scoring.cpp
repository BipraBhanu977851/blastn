#include "scoring.h"
#include <algorithm>
#include <cmath>

// Perform ungapped extension from a seed position
// Match: +2, Mismatch: -1
// Extends both left and right from seed
ExtensionResult extendUngapped(
    const std::string& db_seq,
    const std::string& query,
    int db_seed_pos,
    int q_seed_pos
) {
    ExtensionResult result;
    
    int db_len = static_cast<int>(db_seq.length());
    int q_len = static_cast<int>(query.length());
    
    // Start from seed position (k-mer match, so initial score is k*2)
    // For simplicity, we'll extend and calculate score from scratch
    
    // Extend to the right
    int db_pos = db_seed_pos;
    int q_pos = q_seed_pos;
    int right_score = 0;
    int best_right_score = 0;
    int best_right_end_db = db_seed_pos;
    int best_right_end_q = q_seed_pos;
    
    // Extend right until we can't extend further or score drops too much
    while (db_pos + 1 < db_len && q_pos + 1 < q_len) {
        db_pos++;
        q_pos++;
        
        if (db_seq[db_pos] == query[q_pos]) {
            right_score += 2;  // Match
        } else {
            right_score -= 1;  // Mismatch
        }
        
        // Keep track of best cumulative score
        if (right_score > best_right_score) {
            best_right_score = right_score;
            best_right_end_db = db_pos;
            best_right_end_q = q_pos;
        }
        
        // Stop if score becomes too negative (drop threshold)
        // This prevents extending through very poor regions
        if (right_score < best_right_score - 20) {
            break;
        }
    }
    
    // Extend to the left
    int left_score = 0;
    int best_left_score = 0;
    int best_left_start_db = db_seed_pos;
    int best_left_start_q = q_seed_pos;
    
    db_pos = db_seed_pos;
    q_pos = q_seed_pos;
    
    while (db_pos > 0 && q_pos > 0) {
        db_pos--;
        q_pos--;
        
        if (db_seq[db_pos] == query[q_pos]) {
            left_score += 2;  // Match
        } else {
            left_score -= 1;  // Mismatch
        }
        
        // Keep track of best cumulative score
        if (left_score > best_left_score) {
            best_left_score = left_score;
            best_left_start_db = db_pos;
            best_left_start_q = q_pos;
        }
        
        // Stop if score becomes too negative
        if (left_score < best_left_score - 20) {
            break;
        }
    }
    
    // Combine left and right extensions
    // The total score is left_score + right_score (seed position counted once in each)
    // But since we're extending from seed, we need to account for the seed itself
    // For simplicity, we'll recalculate the full alignment score
    result.db_start = best_left_start_db;
    result.db_end = best_right_end_db;
    result.q_start = best_left_start_q;
    result.q_end = best_right_end_q;
    
    // Recalculate total score for the entire alignment
    int total_score = 0;
    db_pos = result.db_start;
    q_pos = result.q_start;
    while (db_pos <= result.db_end && q_pos <= result.q_end) {
        if (db_seq[db_pos] == query[q_pos]) {
            total_score += 2;
        } else {
            total_score -= 1;
        }
        db_pos++;
        q_pos++;
    }
    result.score = total_score;
    
    // Calculate percent identity
    result.identity = calculateIdentity(
        db_seq, query,
        result.db_start, result.db_end,
        result.q_start, result.q_end
    );
    
    return result;
}

// Calculate percent identity for an alignment
double calculateIdentity(
    const std::string& db_seq,
    const std::string& query,
    int db_start,
    int db_end,
    int q_start,
    int q_end
) {
    if (db_end < db_start || q_end < q_start) {
        return 0.0;
    }
    
    int matches = 0;
    int total = 0;
    
    int db_pos = db_start;
    int q_pos = q_start;
    
    while (db_pos <= db_end && q_pos <= q_end) {
        total++;
        if (db_seq[db_pos] == query[q_pos]) {
            matches++;
        }
        db_pos++;
        q_pos++;
    }
    
    if (total == 0) return 0.0;
    
    return (100.0 * matches) / total;
}

