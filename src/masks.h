#pragma once

#include <string>
#include <iostream>
#include <glpk.h>

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
    of << std::endl;
}

/// For the given masked superstring output the same superstring with mask with minimal number of runs of ones.
void OptimizeRuns(std::string path, kh_S64_t *kMers, std::ostream &of, int k, bool complements) {
    kh_O64_t *intervals = kh_init_O64();
    auto [size, rows] = ReadIntervals(intervals, kMers, path, k, complements, of, nullptr);
    glp_prob *lp;
    lp = glp_create_prob();
    auto *ia = new int[size + 1];
    auto *ja = new int[size + 1];
    auto *ar = new double[size + 1];
    glp_set_obj_dir(lp, GLP_MIN);
    glp_add_rows(lp, rows);
    for (size_t i = 0; i < rows; ++i) {
        glp_set_row_bnds(lp, i + 1, GLP_LO, 1.0, 0.0);
    }
    glp_add_cols(lp, kh_size(intervals));
    for (size_t i = 0; i < kh_size(intervals); ++i) {
        glp_set_col_bnds(lp, i + 1, GLP_DB, 0.0, 1.0);
        glp_set_obj_coef(lp, i + 1, 1.0);
    }
    size_t index = 0;
    size_t kMer = 0;
    for (auto i = kh_begin(intervals); i != kh_end(intervals); ++i) {
        if (!kh_exist(intervals, i)) continue;
        auto key = kh_key(intervals, i);
        for (auto j : kh_value(intervals, key)) {
            ia[index + 1] = j + 1;
            ja[index + 1] = kMer + 1;
            ar[index + 1] = 1.0;
            ++index;
        }
        ++kMer;
    }
    glp_load_matrix(lp, index, ia, ja, ar);
    glp_simplex(lp, nullptr);

    bool *intervalsSet = new bool[rows];

    for (size_t i = 0; i < rows; ++i) {
        intervalsSet[i] = glp_get_col_prim(lp, i + 1) > 0.5;
    }

    ReadIntervals(nullptr, kMers, path, k, complements, of, intervalsSet);
    of << std::endl;
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
        } else if (algorithm == "runs") {
            OptimizeRuns(path, kMers, of, k, complements);
        } else {
            std::cerr << "Algorithm '" + algorithm + "' not recognized." << std::endl;
            in.close();
            return 1;
        }
        in.close();
    } else {
        throw std::invalid_argument("couldn't open file " + path);
    }
    return 0;
}