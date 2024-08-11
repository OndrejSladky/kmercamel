#pragma once

#include <string>
#include <iostream>
#include <glpk.h>

#include "parser.h"
#include "khash_utils.h"
#include "kmers.h"

inline char Masked(char c, bool mask) {
    int masked_difference = (c <= 'Z') - (int) mask;
    return c + (char) masked_difference * ('a' - 'A');
}

/// For the given masked superstring output the same superstring with mask with minimal/maximal number of ones.
void OptimizeOnes(kseq_t masked_superstring, std::ostream &of, kh_S64_t *kMers, int k, bool complements, bool minimize) {
    kmer_t currentKMer = 0, reverseComplement = 0;
    kmer_t mask = (1 << (2 * k)) - 1;
    kmer_t shift = 2 * (k - 1);
    of << masked_superstring.name.s << std::endl;
    uint8_t ms_validation = 0;
    for (size_t i = 0; i < masked_superstring.seq.l; ++i) {
        auto data = nucleotideToInt[(uint8_t) masked_superstring.seq.s[i]];
        ms_validation |= data;
        currentKMer = ((currentKMer << 2) | data) & mask;
        reverseComplement = (reverseComplement >> 2) | ((kmer_t(3 ^ data)) << shift);
        if (i >= k - 1) {
            kmer_t canonical = ((!complements) || currentKMer < reverseComplement) ? currentKMer : reverseComplement;
            auto kmer_pointer = kh_get_S64(kMers, canonical);
            bool contained = kmer_pointer != kh_end(kMers);
            of << Masked(masked_superstring.seq.s[i - k + 1], contained);
            // If minimizing, erase the k-mer once set.
            if (minimize && contained) {
                kh_del(S64, kMers, kmer_pointer);
            }
        }
    }
    // Print the remaining k-1 characters.
    for (size_t i = 0; i < k - 1; ++i) {
        of << Masked(masked_superstring.seq.s[masked_superstring.seq.l - k + 1 + i], false);
    }
    if (ms_validation >= 4) {
        throw std::invalid_argument("Masked superstring contains invalid characters.");
    }
}

/// Pre-solve to the ILP in minimizing the number of runs.
/// Set the intervals with single-occurring k-mers to 1 and remove k-mers from those intervals.
std::pair<std::vector<int>, std::vector<int>> HeuristicPreSolve(std::vector<std::list<size_t>> &intervalsForKMer, size_t intervalCount, int &mappedSize, size_t &totalIntervals, int &newIntervals) {
    std::vector<int> mapping(intervalsForKMer.size());
    std::vector<int> intervalMapping(intervalCount);

    // Set undecided intervals with a single-occurring k-mer to 1s.
    for (size_t i = 0; i < intervalsForKMer.size(); ++i) {
        if (intervalsForKMer[i].size() == 1) {
            mapping[i] = -1;
            intervalMapping[intervalsForKMer[i].front()] = -1;
        }
    }

    // Remove k-mers from decided intervals.
    for (size_t i = 0; i < intervalsForKMer.size(); ++i) {
        for (auto j : intervalsForKMer[i]) {
            if (intervalMapping[j] == -1) {
                mapping[i] = -1;
            }
        }
    }

    totalIntervals = 0;
    int newIndex = 0;
    for (size_t i = 0; i < intervalsForKMer.size(); ++i) {
        if (mapping[i] == -1) continue;
        mapping[i] = newIndex++;
        totalIntervals += intervalsForKMer[i].size();
    }
    mappedSize = newIndex;

    int newIntervalIndex = 0;
    for (size_t i = 0; i < intervalCount; ++i) {
        if (intervalMapping[i] == -1) continue;
        intervalMapping[i] = newIntervalIndex++;
    }
    newIntervals = newIntervalIndex;


    return {mapping, intervalMapping};
}


/// For the given masked superstring output the same superstring with mask with minimal number of runs of ones.
void OptimizeRuns(std::string path, kh_S64_t *kMers, std::ostream &of, int k, bool complements, bool approximate) {
    kh_O64_t *intervals = kh_init_O64();
    std::vector<std::list<size_t>> intervalsForKMer;
    auto [size, rows] = ReadIntervals(intervals, kMers, intervalsForKMer, path, k, complements, of, nullptr);
    int mappedSize, newIntervals; size_t totalIntervals;
    auto [mapping, intervalMapping] = HeuristicPreSolve(intervalsForKMer, rows, mappedSize, totalIntervals, newIntervals);
    glp_prob *lp;
    lp = glp_create_prob();
    if (mappedSize != 0 && !approximate) {
        auto *ia = new int[totalIntervals + 1];
        auto *ja = new int[totalIntervals + 1];
        auto *ar = new double[totalIntervals + 1];
        glp_set_obj_dir(lp, GLP_MIN);
        // Add a row per each k-mer.
        glp_add_rows(lp, mappedSize);
        for (int i = 0; i < mappedSize; ++i) {
            glp_set_row_bnds(lp, i + 1, GLP_LO, 1.0, 0.0);
        }
        // Add a column per each undecided interval.
        glp_add_cols(lp, newIntervals);
        for (int i = 0; i < newIntervals; ++i) {
            glp_set_col_bnds(lp, i + 1, GLP_LO, 0.0, 1.0);
            glp_set_col_kind(lp, i + 1, GLP_IV);
            glp_set_obj_coef(lp, i + 1, 1.0);
        }
        int index = 0;
        for (auto key = kh_begin(intervals); key != kh_end(intervals); ++key) {
            if (!kh_exist(intervals, key)) continue;
            size_t i = kh_value(intervals, key);
            if (mapping[i] == -1) continue;
            for (auto j: intervalsForKMer[i]) {
                ja[index + 1] = intervalMapping[j] + 1;
                ia[index + 1] = mapping[i] + 1;
                ar[index + 1] = 1.0;
                ++index;
            }
        }
        // Supress glpk output.
        glp_term_out(GLP_OFF);

        glp_load_matrix(lp, index, ia, ja, ar);
        glp_simplex(lp, nullptr);
    }

    bool *intervalsSet = new bool[rows];

    for (size_t i = 0; i < rows; ++i) {
        if (intervalMapping[i] == -1) intervalsSet[i] = true;
        else if (approximate) intervalsSet[i] = mappedSize != 0;
        else intervalsSet[i] = mappedSize == 0 ? false : (glp_get_col_prim(lp, intervalMapping[i] + 1) > 0.5);
    }

    ReadIntervals(nullptr, kMers, intervalsForKMer, path, k, complements, of, intervalsSet);
    of << std::endl;
}

int Optimize(std::string &algorithm, std::string path, std::ostream &of,  int k, bool complements) {
    kseq_t masked_superstring = ReadMaskedSuperstring(path);
    kh_S64_t *kMers = kh_init_S64();
    AddKMers(kMers, masked_superstring.seq.l, masked_superstring.seq.s, k, complements, true);

    if (algorithm == "ones") {
        OptimizeOnes(masked_superstring, of, kMers, k, complements, false);
    } else if (algorithm == "zeros") {
        OptimizeOnes(masked_superstring, of, kMers, k, complements, true);
    } else if (algorithm == "runs") {
        OptimizeRuns(path, kMers, of, k, complements, false);
    } else if (algorithm == "runsapprox") {
        OptimizeRuns(path, kMers, of, k, complements, true);
    } else {
        std::cerr << "Algorithm '" + algorithm + "' not recognized." << std::endl;
        kseq_destroy(&masked_superstring);
        return 1;
    }
    kseq_destroy(&masked_superstring);
    return 0;
}