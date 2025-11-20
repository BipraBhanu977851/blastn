#include "fasta.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

// Parse database FASTA file with multiple sequences
std::vector<Sequence> parseDatabase(const std::string& filename) {
    std::vector<Sequence> database;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open database file: " << filename << std::endl;
        return database;
    }
    
    std::string line;
    Sequence current_seq;
    int seq_index = 0;
    
    while (std::getline(file, line)) {
        // Remove carriage return if present
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        
        if (line.empty()) continue;
        
        // Header line starts with '>'
        if (line[0] == '>') {
            // Save previous sequence if exists
            if (!current_seq.seq.empty()) {
                current_seq.index = seq_index++;
                database.push_back(current_seq);
                current_seq = Sequence();
            }
            
            // Parse header: >id|species
            std::string header = line.substr(1); // Remove '>'
            size_t pipe_pos = header.find('|');
            
            if (pipe_pos != std::string::npos) {
                current_seq.id = header.substr(0, pipe_pos);
                current_seq.species = header.substr(pipe_pos + 1);
            } else {
                // If no pipe, use entire header as ID
                current_seq.id = header;
                current_seq.species = "Unknown";
            }
        } else {
            // Sequence line - convert to uppercase and append
            std::string seq_line = line;
            std::transform(seq_line.begin(), seq_line.end(), seq_line.begin(), ::toupper);
            current_seq.seq += seq_line;
        }
    }
    
    // Don't forget the last sequence
    if (!current_seq.seq.empty()) {
        current_seq.index = seq_index++;
        database.push_back(current_seq);
    }
    
    file.close();
    return database;
}

// Parse query FASTA file (single sequence)
std::string parseQuery(const std::string& filename) {
    std::ifstream file(filename);
    std::string query;
    
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open query file: " << filename << std::endl;
        return query;
    }
    
    std::string line;
    bool in_sequence = false;
    
    while (std::getline(file, line)) {
        // Remove carriage return if present
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        
        if (line.empty()) continue;
        
        // Skip header line
        if (line[0] == '>') {
            in_sequence = true;
            continue;
        }
        
        if (in_sequence) {
            // Convert to uppercase and append
            std::string seq_line = line;
            std::transform(seq_line.begin(), seq_line.end(), seq_line.begin(), ::toupper);
            query += seq_line;
        }
    }
    
    file.close();
    return query;
}

// Parse query FASTA file with multiple sequences
std::vector<Query> parseQueries(const std::string& filename) {
    std::vector<Query> queries;
    std::ifstream file(filename);
    
    if (!file.is_open()) {
        std::cerr << "Error: Cannot open query file: " << filename << std::endl;
        return queries;
    }
    
    std::string line;
    Query current_query;
    
    while (std::getline(file, line)) {
        // Remove carriage return if present
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        
        if (line.empty()) continue;
        
        // Header line starts with '>'
        if (line[0] == '>') {
            // Save previous query if exists
            if (!current_query.seq.empty()) {
                queries.push_back(current_query);
                current_query = Query();
            }
            
            // Parse header: >name
            std::string header = line.substr(1); // Remove '>'
            // Remove any pipe and everything after it
            size_t pipe_pos = header.find('|');
            if (pipe_pos != std::string::npos) {
                current_query.name = header.substr(0, pipe_pos);
            } else {
                current_query.name = header.empty() ? "Unknown" : header;
            }
        } else {
            // Sequence line - convert to uppercase and append
            std::string seq_line = line;
            std::transform(seq_line.begin(), seq_line.end(), seq_line.begin(), ::toupper);
            current_query.seq += seq_line;
        }
    }
    
    // Don't forget the last query
    if (!current_query.seq.empty()) {
        queries.push_back(current_query);
    }
    
    file.close();
    return queries;
}

