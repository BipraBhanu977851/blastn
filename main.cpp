#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>
#include <fstream>
#include <sstream>
#include "fasta.h"
#include "index.h"
#include "search.h"

void printUsage(const char* program_name) {
    std::cerr << "Usage: " << program_name 
              << " --db <database.fasta> --query <query.fasta> [--k <kmer_size>] [--top <N>]"
              << std::endl;
    std::cerr << "Options:" << std::endl;
    std::cerr << "  --db    : Database FASTA file (required)" << std::endl;
    std::cerr << "  --query : Query FASTA file (required)" << std::endl;
    std::cerr << "  --k     : K-mer size (default: 11)" << std::endl;
    std::cerr << "  --top   : Number of top hits per query (default: 2, 0 = all)" << std::endl;
}

// Format range string
std::string formatRange(int start, int end) {
    return std::to_string(start) + "-" + std::to_string(end);
}

// Wrap alignment lines to max 80 characters
void printWrappedAlignment(const std::string& db_seq, const std::string& match_line, 
                          const std::string& q_seq) {
    const int MAX_LINE = 80;
    const int PREFIX_LEN = 6; // "DB:   " or "      " or "Q:   "
    
    int chunk_size = MAX_LINE - PREFIX_LEN;
    size_t len = db_seq.length();
    
    for (size_t i = 0; i < len; i += chunk_size) {
        size_t end = std::min(i + chunk_size, len);
        
        std::cout << "DB:   " << db_seq.substr(i, end - i) << std::endl;
        std::cout << "      " << match_line.substr(i, end - i) << std::endl;
        std::cout << "Q:    " << q_seq.substr(i, end - i) << std::endl;
        
        if (end < len) {
            std::cout << std::endl;
        }
    }
}

int main(int argc, char* argv[]) {
    std::string db_file;
    std::string query_file;
    int k = 11;
    int top_n = 2;  // Default to showing top 2 hits (0 = all)
    
    // Parse command-line arguments
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        
        if (arg == "--db" && i + 1 < argc) {
            db_file = argv[++i];
        } else if (arg == "--query" && i + 1 < argc) {
            query_file = argv[++i];
        } else if (arg == "--k" && i + 1 < argc) {
            k = std::stoi(argv[++i]);
            if (k < 1 || k > 16) {
                std::cerr << "Error: k must be between 1 and 16" << std::endl;
                return 1;
            }
        } else if (arg == "--top" && i + 1 < argc) {
            top_n = std::stoi(argv[++i]);
            if (top_n < 0) {
                std::cerr << "Error: top must be non-negative" << std::endl;
                return 1;
            }
        } else if (arg == "--help" || arg == "-h") {
            printUsage(argv[0]);
            return 0;
        }
    }
    
    // Check required arguments
    if (db_file.empty() || query_file.empty()) {
        std::cerr << "Error: --db and --query are required" << std::endl;
        printUsage(argv[0]);
        return 1;
    }
    
    // Step 1: Parse FASTA files
    std::vector<Sequence> database = parseDatabase(db_file);
    if (database.empty()) {
        std::cerr << "Error: No sequences found in database file" << std::endl;
        return 1;
    }
    
    std::vector<Query> queries = parseQueries(query_file);
    if (queries.empty()) {
        std::cerr << "Error: No queries found in query file" << std::endl;
        return 1;
    }
    
    // Step 2: Build k-mer index
    KmerIndex index = buildIndex(database, k);
    
    // Process each query
    for (size_t q_idx = 0; q_idx < queries.size(); ++q_idx) {
        const Query& query = queries[q_idx];
        
        if (query.seq.empty()) {
            std::cerr << "Warning: Query " << query.name << " is empty, skipping" << std::endl;
            continue;
        }
        
        // Step 3: Search for HSPs
        std::vector<HSP> hsps = findHSPs(query.seq, database, index, k);
        
        // Step 4: Merge overlapping HSPs
        std::vector<HSP> merged_hsps = mergeHSPs(hsps);
        
        // Step 5: Sort by score (descending), then by identity (descending)
        std::sort(merged_hsps.begin(), merged_hsps.end(),
            [](const HSP& a, const HSP& b) {
                if (a.score != b.score) {
                    return a.score > b.score;
                }
                return a.identity > b.identity;
            });
        
        // Step 6: Display results in compact format
        if (merged_hsps.empty()) {
            std::cout << "QUERY: " << query.name << "   (" << query.seq.length()
                      << " bp)" << std::endl;
            std::cout << std::endl;
            std::cout << "BEST HIT: No hits found" << std::endl;
            if (q_idx < queries.size() - 1) {
                std::cout << std::endl;
            }
            continue;
        }
        
        int display_count = (top_n == 0)
            ? static_cast<int>(merged_hsps.size())
            : std::min(top_n, static_cast<int>(merged_hsps.size()));
        const HSP& best_hsp = merged_hsps[0];
        const Sequence& best_seq = database[best_hsp.sid];
        
        // Print headers
        std::cout << "QUERY: " << query.name << "   (" << query.seq.length()
                  << " bp)" << std::endl;
        std::cout << std::endl;
        std::cout << "BEST HIT: " << best_seq.species << std::endl;
        std::cout << std::endl;
        
        // Print summary table
        const std::string table_header =
            "Species        Score   Identity   DB Range   Q Range";
        std::cout << table_header << std::endl;
        std::cout << std::string(table_header.length(), '-') << std::endl;
        
        for (int i = 0; i < display_count; ++i) {
            const HSP& hsp = merged_hsps[i];
            const Sequence& seq = database[hsp.sid];
            
            // Truncate species name if too long
            std::string species_display = seq.species;
            if (species_display.length() > 14) {
                species_display = species_display.substr(0, 11) + "...";
            }
            
            std::ostringstream identity_stream;
            identity_stream << std::fixed << std::setprecision(2) << hsp.identity
                            << "%";
            std::string identity_str = identity_stream.str();
            std::string db_range = formatRange(hsp.db_start, hsp.db_end);
            std::string q_range = formatRange(hsp.q_start, hsp.q_end);
            
            std::cout << std::left << std::setw(14) << species_display;
            std::cout << std::right << std::setw(7) << hsp.score;
            std::cout << std::right << std::setw(12) << identity_str;
            std::cout << std::right << std::setw(11) << db_range;
            std::cout << std::right << std::setw(9) << q_range << std::endl;
            std::cout << std::left;
        }
        
        std::cout << std::endl;
        
        // Print alignment blocks
        for (int i = 0; i < display_count; ++i) {
            const HSP& hsp = merged_hsps[i];
            const Sequence& seq = database[hsp.sid];
            
            if (display_count > 1) {
                std::cout << "Hit #" << (i + 1) << " (" << seq.species << ")"
                          << std::endl;
            }
            
            std::string alignment = getAlignment(
                seq.seq, query.seq,
                hsp.db_start, hsp.db_end,
                hsp.q_start, hsp.q_end
            );
            
            if (!alignment.empty()) {
                size_t first_nl = alignment.find('\n');
                size_t second_nl = alignment.find('\n', first_nl + 1);
                if (first_nl != std::string::npos &&
                    second_nl != std::string::npos) {
                    std::string db_seq = alignment.substr(0, first_nl);
                    std::string match_line =
                        alignment.substr(first_nl + 1,
                                         second_nl - first_nl - 1);
                    std::string q_seq = alignment.substr(second_nl + 1);
                    
                    printWrappedAlignment(db_seq, match_line, q_seq);
                }
            }
            
            if (i < display_count - 1) {
                std::cout << std::endl;
            }
        }
        
        // Add separator between queries
        if (q_idx < queries.size() - 1) {
            std::cout << std::endl;
        }
    }
    
    return 0;
}
