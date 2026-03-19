#pragma once

#include "global.h"
#include "global_sparse.h"
#include "kmers.h"

/// Return the length of the cycle cover which lower bounds the superstring length.
template <typename kmer_t, typename kh_wrapper_t>
size_t LowerBoundLength(kh_wrapper_t wrapper, kmer_t kmer_type, std::vector<simplitig_t> simplitigs, int k, bool complements) {
    auto cycle_cover = OverlapHamiltonianPath(wrapper, kmer_type, simplitigs, k, complements, true);
    WriteLog("Finished 2. part: Hamiltonian path.");
    size_t res = 0;
    for (auto &simplitig : simplitigs) {
        res += simplitig.size() / (2 - complements);
    }
    for (auto &overlap : cycle_cover.second) {
        res -= overlap;
    }
    res /= (1 + complements);
    WriteLog("Finished 3. part: lower bound = " + std::to_string(res) + ".");
    return res;
}

/// Same as LowerBoundLength for the k-mer overlap graph (as used by GlobalSparse / PartialPreSort).
template <typename kmer_t, typename kh_wrapper_t>
size_t LowerBoundLengthSparse(kh_wrapper_t wrapper, std::vector<kmer_t> &kMers, int k, bool complements) {
    auto cycle_cover = OverlapHamiltonianPathSparse(wrapper, kMers, k, complements, true);
    WriteLog("Finished 2. part: Hamiltonian path.");
    size_t n = kMers.size();
    size_t res = n * k * (1 + complements);
    for (auto &overlap : cycle_cover.second) {
        res -= overlap;
    }
    res /= (1 + complements);
    WriteLog("Finished 3. part: lower bound = " + std::to_string(res) + ".");
    return res;
}