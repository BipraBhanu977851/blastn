#include "search.h"
#include "index.h"
#include <algorithm>
#include <set>

// Find all HSPs for a query sequence
std::vector<HSP> findHSPs(
    const std::string& query,
    const std::vector<Sequence>& database,
    const KmerIndex& index,
    int k
) {
    std::vector<HSP> hsps;
    int q_len = static_cast<int>(query.length());
    
    // For each k-mer in query
    for (int q_pos = 0; q_pos <= q_len - k; ++q_pos) {
        // Get k-mer at query position
        uint32_t kmer_key = getKmerAt(query, q_pos, k);
        
        // Skip invalid k-mers
        if (kmer_key == 0) {
            bool valid = true;
            for (int j = 0; j < k; ++j) {
                char c = query[q_pos + j];
                if (c != 'A' && c != 'C' && c != 'G' && c != 'T' &&
                    c != 'a' && c != 'c' && c != 'g' && c != 't') {
                    valid = false;
                    break;
                }
            }
            if (!valid) continue;
        }
        
        // Lookup in hash table
        auto it = index.find(kmer_key);
        if (it != index.end()) {
            // For each hit in database
            for (const auto& hit : it->second) {
                int db_seq_idx = hit.first;
                int db_seed_pos = hit.second;
                
                // Perform ungapped extension
                ExtensionResult ext = extendUngapped(
                    database[db_seq_idx].seq,
                    query,
                    db_seed_pos,
                    q_pos
                );
                
                // Create HSP
                HSP hsp;
                hsp.sid = db_seq_idx;
                hsp.db_start = ext.db_start;
                hsp.db_end = ext.db_end;
                hsp.q_start = ext.q_start;
                hsp.q_end = ext.q_end;
                hsp.score = ext.score;
                hsp.identity = ext.identity;
                
                hsps.push_back(hsp);
            }
        }
    }
    
    return hsps;
}

// Merge overlapping HSPs for the same sequence
// Keeps the best scoring HSP when overlaps occur
std::vector<HSP> mergeHSPs(const std::vector<HSP>& hsps) {
    if (hsps.empty()) return hsps;
    
    // Group HSPs by sequence ID
    std::vector<std::vector<HSP>> by_sequence;
    std::vector<int> seq_ids;
    
    for (const auto& hsp : hsps) {
        auto it = std::find(seq_ids.begin(), seq_ids.end(), hsp.sid);
        if (it == seq_ids.end()) {
            seq_ids.push_back(hsp.sid);
            by_sequence.push_back({hsp});
        } else {
            int idx = std::distance(seq_ids.begin(), it);
            by_sequence[idx].push_back(hsp);
        }
    }
    
    std::vector<HSP> merged;
    
    // For each sequence, merge overlapping HSPs
    for (auto& seq_hsps : by_sequence) {
        // Sort by database start position
        std::sort(seq_hsps.begin(), seq_hsps.end(),
            [](const HSP& a, const HSP& b) {
                return a.db_start < b.db_start;
            });
        
        // Simple merging: keep non-overlapping or best scoring overlapping
        for (size_t i = 0; i < seq_hsps.size(); ++i) {
            bool overlap = false;
            
            // Check if this HSP overlaps with any already merged for this sequence
            for (const auto& m : merged) {
                if (m.sid == seq_hsps[i].sid) {
                    // Check overlap: ranges overlap if not (end1 < start2 || end2 < start1)
                    if (!(seq_hsps[i].db_end < m.db_start || m.db_end < seq_hsps[i].db_start)) {
                        overlap = true;
                        // If current HSP has better score, we'll replace (handled below)
                        break;
                    }
                }
            }
            
            if (!overlap) {
                merged.push_back(seq_hsps[i]);
            } else {
                // Replace overlapping HSP if this one is better
                for (auto& m : merged) {
                    if (m.sid == seq_hsps[i].sid) {
                        if (!(seq_hsps[i].db_end < m.db_start || m.db_end < seq_hsps[i].db_start)) {
                            if (seq_hsps[i].score > m.score ||
                                (seq_hsps[i].score == m.score && seq_hsps[i].identity > m.identity)) {
                                m = seq_hsps[i];
                            }
                            break;
                        }
                    }
                }
            }
        }
    }
    
    return merged;
}

// Get alignment string representation
std::string getAlignment(
    const std::string& db_seq,
    const std::string& query,
    int db_start,
    int db_end,
    int q_start,
    int q_end
) {
    std::string alignment;
    
    if (db_end < db_start || q_end < q_start) {
        return alignment;
    }
    
    int db_pos = db_start;
    int q_pos = q_start;
    
    std::string db_line;
    std::string match_line;
    std::string q_line;
    
    while (db_pos <= db_end && q_pos <= q_end) {
        char db_char = db_seq[db_pos];
        char q_char = query[q_pos];
        
        db_line += db_char;
        q_line += q_char;
        
        if (db_char == q_char) {
            match_line += '|';
        } else {
            match_line += ' ';
        }
        
        db_pos++;
        q_pos++;
    }
    
    alignment = db_line + "\n" + match_line + "\n" + q_line;
    
    return alignment;
}

