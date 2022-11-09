#pragma once
#include <string>
#include <unordered_set>
#include <fstream>
#include <algorithm>

#include "models.h"

/// Record of one fasta sequence.
struct FastaRecord {
    // The name of the record. Starts with '>'.
    std::string name;
    // The genetic sequence.
    std::string sequence;
};

/// Read fasta file with given path.
std::vector<FastaRecord> ReadFasta(std::string &path) {
    std::ifstream fasta(path);
    std::vector<FastaRecord> records;
    std::string line;
    if (fasta.is_open()) {
        while (std::getline(fasta, line)) {
            // Start new record.
            if (!line.empty() && line[0] == '>') {
                records.push_back({line, ""});
            // Ignore lines before the first record.
            } else if (!records.empty()) {
                // Append to the last record.
                records[records.size() - 1].sequence += line;
            }
        }
        fasta.close();
    } else {
        throw std::invalid_argument("couldn't open file " + path);
    }
    return records;
}

/// Create a list of unique k-mers in no particular order.
/// This runs in O(k*data.size) expected time.
void AddKMersFromSequence(std::unordered_set<std::string> &kMers, std::string data, int k) {
    // Convert the sequence to uppercase letters.
    std::transform(data.begin(), data.end(), data.begin(), toupper);
    size_t possibleKMerEnd = k;
    for (size_t i = k; i <= data.size(); ++i) {
        if (data[i-1] != 'A' && data[i-1] != 'C' && data[i-1] != 'G' && data[i-1] != 'T') {
            // Skip this and the next k-1 k-mers.
            possibleKMerEnd = i + k;
        }
        if (i >= possibleKMerEnd) {
            kMers.insert(data.substr(i - k, k));
        }
    }
}

/// Create a list of unique k-mers in no particular order.
/// This runs in O(k*data.size) expected time.
std::vector<KMer> ConstructKMers(std::vector<FastaRecord> &data, int k) {
    std::unordered_set<std::string> uniqueKMers;
    for (auto &&record : data) {
        AddKMersFromSequence(uniqueKMers, record.sequence, k);
    }
    std::vector<KMer> result;
    for (auto &kMer : uniqueKMers) {
        result.push_back(KMer{kMer});
    }
    return result;
}
