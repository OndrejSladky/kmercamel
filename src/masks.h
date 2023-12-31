#pragma once

#include <string>
#include <iostream>

#include "parser.h"
#include "khash_utils.h"
#include "kmers.h"

/// For the given masked superstring output the same superstring with mask with minimal/maximal number of ones.
void OptimizeOnes(std::ifstream &in, std::ostream &of, kh_S64_t *kMers, int k, bool complements, bool minimize) {
    char c;
    int beforeKMerEnd = k;
    kmer_t currentKMer = 0;
    kmer_t mask = (((kmer_t) 1) <<  (2 * k) ) - 1;
    bool readingHeader = false;
    while (in >> std::noskipws >> c) {
        // Start of a header.
        if (c == '>') {
            readingHeader = true;
            currentKMer &= mask >> 1;
            for (char ch : NumberToKMer(currentKMer, k - beforeKMerEnd)) of << UpperToLower(ch);
            currentKMer = 0;
            beforeKMerEnd = k;
        }
        // Reprint the header.
        if (readingHeader) of << c;
        if (c == '\n') readingHeader = false;
        if (readingHeader) continue;
        auto data = NucleotideToInt(c);
        // Disregard white space.
        if (c == '\n' || c == '\r' || c == ' ') continue;
        // Break at an unknown character.
        if (data == -1) {
            currentKMer &= mask >> 1;
            for (char ch : NumberToKMer(currentKMer, k - beforeKMerEnd)) of << UpperToLower(ch);
            currentKMer = 0;
            beforeKMerEnd = k;
            continue;
        }
        currentKMer <<= 2;
        currentKMer &= mask;
        currentKMer |= data;
        --beforeKMerEnd;
        kmer_t reverseComplement = ReverseComplement(currentKMer, k);
        if (beforeKMerEnd == 0) {
            ++beforeKMerEnd;
            char to_print = NucleotideAtIndex(currentKMer, k, 0);
            auto key = kh_get_S64(kMers, currentKMer);
            if (key != kh_end(kMers)) {
                if (minimize) {
                    kh_del_S64(kMers, key);
                }
            } else if (complements && (key = kh_get_S64(kMers, reverseComplement)) != kh_end(kMers)) {
                if (minimize) {
                    kh_del_S64(kMers, key);
                }
            } else {
                to_print = UpperToLower(to_print);
            }
            of << to_print;
        }
    }
    currentKMer &= mask >> 1;
    for (char ch : NumberToKMer(currentKMer, k - beforeKMerEnd)) of << UpperToLower(ch);
    of << "\n";
}

int Optimize(std::string &algorithm, std::string path, std::ostream &of,  int k, bool complements) {
    kh_S64_t *kMers = kh_init_S64();
    ReadKMers(kMers, path, k, complements, true);

    std::ifstream in(path);
    if (in.is_open()) {
        if (algorithm == "ones") {
            OptimizeOnes(in, of, kMers, k, complements, false);
        } else if (algorithm == "zeros") {
            OptimizeOnes(in, of, kMers, k, complements, true);
        } else {
            if (algorithm == "runs") {
                std::cerr
                        << "Minimization of the number of runs of ones is not supported. Please use the Python script from supplement."
                        << std::endl;
            } else {
                std::cerr << "Algorithm '" + algorithm + "' not recognized." << std::endl;
            }
            in.close();
            return 1;
        }
        in.close();
    } else {
        throw std::invalid_argument("couldn't open file " + path);
    }
    return 0;
}