#pragma once

#include "kmers.h"
#include "khash.h"

#include <vector>

KHASH_SET_INIT_INT64(S64)

/// Determine whether the k-mer or its reverse complement is present.
bool containsKMer(kh_S64_t *kMers, kmer_t kMer, int k, bool complements) {
    bool ret = kh_get_S64(kMers, kMer) != kh_end(kMers);
    if (complements) ret |= kh_get_S64(kMers, ReverseComplement(kMer, k )) != kh_end(kMers);
    return ret;
}

/// Remove the k-mer and its reverse complement.
void eraseKMer(kh_S64_t *kMers, kmer_t kMer, int k, bool complements) {
    auto key = kh_get_S64(kMers, kMer);
    if (key != kh_end(kMers)) {
        kh_del_S64(kMers, key);
    }
    if (complements) {
        kmer_t reverseComplement = ReverseComplement(kMer, k);
        key = kh_get_S64(kMers, reverseComplement);
        if (key != kh_end(kMers)) kh_del_S64(kMers, key);
    }
}

/// Return the next k-mer in the k-mer set and update the index.
kmer_t nextKMer(kh_S64_t *kMers, size_t &lastIndex) {
    for (size_t i = kh_begin(kMers) + lastIndex; i != kh_end(kMers); ++i, ++lastIndex) {
        if (!kh_exist(kMers, i)) continue;
        return kh_key(kMers, i);
    }
    // No more k-mers.
    return -1;
}

/// Construct a vector of the k-mer set in an arbitrary order.
std::vector<kmer_t> kMersToVec(kh_S64_t *kMers) {
    std::vector<kmer_t> res(kh_size(kMers));
    size_t index = 0;
    for (auto i = kh_begin(kMers); i != kh_end(kMers); ++i) {
        if (!kh_exist(kMers, i)) continue;
        res[index++] = kh_key(kMers, i);
    }
    return res;
}




