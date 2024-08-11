#pragma once
#include <string>
#include <unordered_set>
#include <fstream>
#include <algorithm>

#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "models.h"
#include "kmers.h"
#include "khash_utils.h"

/// Record of one fasta sequence.
struct FastaRecord {
    // The name of the record. Starts with '>'.
    std::string name;
    // The genetic sequence.
    std::string sequence;
};

/// Read fasta file with given path.
std::vector<FastaRecord> ReadFasta(std::string &path) {
    gzFile fp;
    kseq_t *seq;
    std::vector<FastaRecord> records;

    fp = gzopen(path.c_str(), "r");
    if (fp == nullptr) {
        throw std::invalid_argument("couldn't open file " + path);
    }
    seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        records.push_back(FastaRecord{
            seq->name.s,
            seq->seq.s,
        });
    }
    kseq_destroy(seq);
    gzclose(fp);
    return records;
}

/// Create a list of unique k-mers in no particular order.
/// This runs in O(k*data.size) expected time.
void AddKMersFromSequence(std::unordered_set<std::string> &kMers, std::string &data, int k) {
    // Convert the sequence to uppercase letters.
    std::transform(data.begin(), data.end(), data.begin(), toupper);
    size_t possibleKMerEnd = k;
    for (size_t i = 1; i <= data.size(); ++i) {
        if (data[i-1] != 'A' && data[i-1] != 'C' && data[i-1] != 'G' && data[i-1] != 'T') {
            // Skip this and the next k-1 k-mers.
            possibleKMerEnd = i + k;
        }
        if (i >= possibleKMerEnd) {
            kMers.insert(data.substr(i - k, k));
        }
    }
}

/// Return only the subset of k-mers where no two k-mers are complements of one another.
/// The k-mers are chosen with no particular property.
std::unordered_set<std::string> FilterKMersWithComplement(std::unordered_set<std::string> &kMers) {
    std::unordered_set<std::string> ret;
    for (auto &&x : kMers) {
        auto kMer = KMer{x};
        if (ret.count(ReverseComplement(kMer).value) == 0) ret.insert(x);
    }
    return ret;
}

/// Create a list of unique k-mers in no particular order.
/// This runs in O(k*data.size) expected time.
std::vector<KMer> ConstructKMers(std::vector<FastaRecord> &data, int k, bool complements) {
    std::unordered_set<std::string> uniqueKMers;
    for (auto &&record : data) {
        AddKMersFromSequence(uniqueKMers, record.sequence, k);
    }
    if (complements) uniqueKMers = FilterKMersWithComplement(uniqueKMers);
    std::vector<KMer> result;
    for (auto &kMer : uniqueKMers) {
        result.push_back(KMer{kMer});
    }
    return result;
}

/// Fill the k-mer dictionary with k-mers from the given sequence.
/// If complements is true, always add the canonical k-mers.
/// If case_sensitive is true, add the k-mer only if it starts with an upper case letter.
void AddKMers(kh_S64_t *kMers, size_t sequence_length, const char* sequence, int64_t k, bool complements, bool case_sensitive = false) {
    int64_t currentLength = 0;
    kmer_t currentKMer = 0, reverseComplement = 0;
    kmer_t cases = 0;
    kmer_t mask = (((kmer_t) 1) <<  (2 * k) ) - 1;
    kmer_t shift = 2 * (k - 1);
    for (size_t i = 0; i < sequence_length; ++i) {
        auto data = nucleotideToInt[(uint8_t)sequence[i]];
        if (data >= 4) {
            // Restart if "N"-like nucleotide.
            currentKMer = reverseComplement = 0;
            currentLength = 0;
            continue;
        }
        currentKMer = ((currentKMer << 2) | data) & mask; 
        reverseComplement = (reverseComplement >> 2) | ((kmer_t(3 ^ data)) << shift);
        // K-mer is present if it is upper case or case-insensitive.
        cases = (cases | (!case_sensitive || sequence[i] <= 'Z')) << 1;
        if ((++currentLength >= k) && (cases & (kmer_t(1) << k))) {
            // Add the canonical k-mer to the dictionary.
            kmer_t canonical = ((!complements) || currentKMer < reverseComplement) ? currentKMer : reverseComplement;
            int ret;
            kh_put_S64(kMers, canonical, &ret);
        }
    }
}


/// Load a dictionary of k-mers from a fasta file.
/// If complements is true, add the canonical k-mers.
void ReadKMers(kh_S64_t *kMers, std::string &path, int k, bool complements, bool case_sensitive = false) {
    gzFile fp;
    kseq_t *seq;

    fp = gzopen(path.c_str(), "r");
    if (fp == nullptr) {
        throw std::invalid_argument("couldn't open file " + path);
    }

    seq = kseq_init(fp);
    while (kseq_read(seq) >= 0) {
        AddKMers(kMers, seq->seq.l, seq->seq.s, k, complements, case_sensitive);
    }

    kseq_destroy(seq);
    gzclose(fp);
}

/// Print the k-mer tail that has [beforeKMerEnd] steps to become a full k-mer.
void PrintRemainingKMer(kmer_t currentKMer, int beforeKMerEnd, int k, std::ostream &of) {
    currentKMer <<= 2 * beforeKMerEnd;
    for (int i = 0; i < k - beforeKMerEnd; ++i) {
        char c = NucleotideAtIndex(currentKMer, k, i) - 'A' + 'a';
        of << c;
    }
}

/// Read or set the intervals.
/// If [setIntervals] is provided reprint the given files with the corresponding intervals set to 1.
/// Otherwise, read the intervals in which each k-mer occurs.
std::pair<size_t, size_t> ReadIntervals(kh_O64_t *intervals, kh_S64_t *kMers, std::vector<std::list<size_t>> &intervalsForKmer,
                                        std::string &path, int k, bool complements, std::ostream &of, const bool* setIntervals = nullptr) {
    std::ifstream fasta(path);
    if (fasta.is_open()) {
        bool reading = setIntervals == nullptr;
        size_t occurrences = 0;
        char c;
        int beforeKMerEnd = k;
        kmer_t currentKMer = 0;
        kmer_t mask = (((kmer_t) 1) <<  (2 * k) ) - 1;
        size_t currentInterval = 0;
        bool interval_used = false;
        bool readingHeader = false;
        while (fasta >> std::noskipws >> c) {
            if (c == '>') {
                if (!reading) PrintRemainingKMer(currentKMer, beforeKMerEnd, k, of);
                readingHeader = true;
                currentKMer = 0;
                beforeKMerEnd = k;
                currentInterval += interval_used;
                interval_used = false;
            }
            // Reprint the header.
            if (readingHeader && !reading) of << c;
            if (c == '\n') readingHeader = false;
            if (readingHeader) continue;
            auto data = NucleotideToInt(c);
            // Disregard white space.
            if (c == '\n' || c == '\r' || c == ' ') continue;
            if (data == -1) {
                if (!reading) PrintRemainingKMer(currentKMer, beforeKMerEnd, k, of);
                currentKMer = 0;
                beforeKMerEnd = k;
                currentInterval += interval_used;
                interval_used = false;
                continue;
            }
            currentKMer <<= 2;
            currentKMer &= mask;
            currentKMer |= data;
            --beforeKMerEnd;
            if (beforeKMerEnd == 0) {
                bool represented = containsKMer(kMers, currentKMer, k, complements);
                bool set = false;
                if (represented) {
                    interval_used = true;
                    if (reading) occurrences += appendInterval(intervals, intervalsForKmer, currentKMer, currentInterval, k, complements);
                    else set = setIntervals[currentInterval];
                } else {
                    currentInterval += interval_used;
                    interval_used = false;
                }
                if (!reading) {
                    char toPrint = NucleotideAtIndex(currentKMer, k, 0);
                    if (set && toPrint >= 'a') {
                        toPrint += 'A' - 'a';
                    } else if (!set && toPrint <= 'Z') {
                        toPrint -= 'A' - 'a';
                    }
                    of << toPrint;
                }
                beforeKMerEnd++;
            }
        }
        if (!reading) PrintRemainingKMer(currentKMer, beforeKMerEnd, k, of);
        fasta.close();
        return {occurrences, currentInterval + interval_used};
    } else {
        throw std::invalid_argument("couldn't open file " + path);
    }
}
