# Simple BLASTN Program

A simple BLASTN-like program in C++ that identifies unknown species from a query DNA sequence using seed-and-extend algorithm with k-mer indexing.

## Overview

This program implements a simplified version of BLASTN (Basic Local Alignment Search Tool for Nucleotides) that:

- Reads a database FASTA file with multiple sequences labeled as `>id|species`
- Reads a query FASTA file (single sequence)
- Builds a k-mer hash index from the database (default k=11)
- Uses ungapped local extension for scoring (+2 match, -1 mismatch)
- Outputs the top N species with highest alignment score and percent identity

## Algorithm Explanation

### 1. K-mer Indexing

**K-mer Hashing:**
- Each k-mer (substring of length k) is encoded using 2-bit encoding:
  - A → 00 (0)
  - C → 01 (1)
  - G → 10 (2)
  - T → 11 (3)
- This allows efficient integer-based hashing (supports k up to 16 with 32-bit integers)
- Example: k-mer "ACGT" = 0×2⁶ + 1×2⁴ + 2×2² + 3×2⁰ = 0 + 16 + 8 + 3 = 27

**Rolling Hash:**
- The program uses a sliding window approach to extract all k-mers from database sequences
- For each position i in a sequence, the k-mer starting at position i is encoded and stored
- This creates a hash table mapping k-mer keys to lists of (sequence_index, position) pairs

**Hash Table Structure:**
```cpp
unordered_map<uint32_t, vector<pair<int, int>>>
```
- Key: Encoded k-mer (uint32_t)
- Value: List of (database_sequence_index, position) pairs where this k-mer appears

### 2. Seed Extension

**Seed Finding:**
- For each k-mer in the query sequence, the program looks it up in the hash table
- Each match is a "seed" - a potential starting point for alignment

**Ungapped Extension:**
- From each seed position, the program extends both left and right without gaps
- Scoring system:
  - Match: +2 points
  - Mismatch: -1 point
- Extension continues until the score drops significantly (threshold-based stopping)
- The extension with the best score is kept as a High Scoring Pair (HSP)

### 3. HSP Management

**HSP Structure:**
```cpp
struct HSP {
    int sid;           // Sequence index in database
    int db_start;      // Start position in database
    int db_end;        // End position in database
    int q_start;       // Start position in query
    int q_end;         // End position in query
    int score;         // Alignment score
    double identity;   // Percent identity
};
```

**Merging Overlapping HSPs:**
- HSPs from the same database sequence that overlap are merged
- When overlaps occur, the HSP with the best score is kept
- This prevents duplicate reporting of the same alignment region

### 4. Ranking and Output

Results are sorted by:
1. Alignment score (descending)
2. Percent identity (descending)

The top N results are displayed with:
- Sequence ID and species name
- Database and query positions
- Alignment score and percent identity
- Visual alignment representation

## Compilation

### Prerequisites
- C++17 compatible compiler (g++ or clang++)
- Make utility

### Build Instructions

```bash
# Compile the program
make

# Or manually:
g++ -std=c++17 -O2 main.cpp fasta.cpp index.cpp search.cpp scoring.cpp -o simple_blastn

# Clean build artifacts
make clean

# Rebuild from scratch
make rebuild
```

## Usage

```bash
./simple_blastn --db <database.fasta> --query <query.fasta> [--k <kmer_size>] [--top <N>]
```

### Arguments

- `--db <file>`: Database FASTA file (required)
  - Format: `>id|species\nsequence`
  - Can contain multiple sequences

- `--query <file>`: Query FASTA file (required)
  - Format: `>query_id\nsequence`
  - Should contain a single sequence

- `--k <size>`: K-mer size (optional, default: 11)
  - Must be between 1 and 16
  - Larger k = fewer false positives, but may miss some matches

- `--top <N>`: Number of top hits to display (optional, default: 5)

### Example

```bash
./simple_blastn --db database.fasta --query query.fasta --k 11 --top 5
```

## Input Format

### Database FASTA (`database.fasta`)
```
>seq1|Escherichia_coli
ATGCTAGCTAGCTTGACCTGATGCTAGCTAGCTAGCTGACTGATCG
>seq2|Bacillus_subtilis
GCTAGCTTGACCGTAGCTAGCTAAAACCCGGGTTTACGATCGATC
>seq3|Saccharomyces_cerevisiae
TTAACCGGTTAGCTAGGCTAGCTAGCTTTGGGCCCATGCTAGCTAG
```

### Query FASTA (`query.fasta`)
```
>query
GCTAGCTTGACCGTAGCTAGCT
```

## Output Format

For each top hit, the program displays:

```
Result #1:
  Sequence ID: seq2
  Species: Bacillus_subtilis
  DB Position: 0-21 (length: 22)
  Query Position: 0-21 (length: 22)
  Alignment Score: 44
  Percent Identity: 100.00%
  Alignment:
    DB:  GCTAGCTTGACCGTAGCTAGCT
         ||||||||||||||||||||
    Q:   GCTAGCTTGACCGTAGCTAGCT
```

## Project Structure

```
.
├── main.cpp          # Main program with command-line interface
├── fasta.h/cpp       # FASTA file parsing functions
├── index.h/cpp       # K-mer indexing and hash table building
├── search.h/cpp      # HSP finding and merging
├── scoring.h/cpp     # Ungapped extension and scoring
├── Makefile          # Build configuration
└── README.md         # This file
```

## Limitations

- **No gap penalties**: Only ungapped alignments are computed
- **No E-values**: Statistical significance is not calculated
- **Simple extension**: Extension stops when score drops, not using dynamic programming
- **Limited k-mer size**: Maximum k=16 due to 32-bit encoding
- **No reverse complement**: Only searches forward strand

## Future Enhancements

Possible improvements:
- Gapped alignment using dynamic programming
- E-value calculation for statistical significance
- Reverse complement search
- Support for larger k-mer sizes
- Parallel processing for large databases
- More sophisticated HSP merging strategies

## License

This is an educational implementation of BLASTN concepts.

