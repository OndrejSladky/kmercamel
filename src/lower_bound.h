#pragma once

#include "global.h"
#include "kmers.h"

/// Return the length of the cycle cover which lower bounds the superstring length.
template <typename kmer_t, typename kh_wrapper_t>
size_t LowerBoundLength(kh_wrapper_t wrapper, kmer_t kmer_type, std::vector<simplitig_t> simplitigs, int k, bool complements) {
    auto cycle_cover = OverlapHamiltonianPath(wrapper, kmer_type, simplitigs, k, complements, true);
    size_t res = 0;
    for (auto &simplitig : simplitigs) {
        res += simplitig.size();
    }
    for (auto &overlap : cycle_cover.second) {
        res -= size_t(overlap);
    }
    return res / (1 + complements);
}