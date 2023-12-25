#pragma once

#include "kmers.h"
#include "khash.h"

#include <unordered_set>

KHASH_MAP_INIT_INT64(P64, size_t)

/// Determine whether the k-mer or its reverse complement is present.
bool containsKMer(kh_P64_t *kMers, int64_t kMer, int k, bool complements) {
    bool ret = kh_get_P64(kMers, kMer) != kh_end(kMers);
    if (complements) ret |= kh_get_P64(kMers, ReverseComplement(kMer, k )) != kh_end(kMers);
    return ret;
}

/// Remove the k-mer and its reverse complement.
void eraseKMer(kh_P64_t *kMers, int64_t kMer, int k, bool complements) {
    auto key = kh_get_P64(kMers, kMer);
    if (key != kh_end(kMers)) kh_del_P64(kMers, key);
    if (complements) {
        int64_t reverseComplement = ReverseComplement(kMer, k);
        key = kh_get_P64(kMers, reverseComplement);
        if (key != kh_end(kMers)) kh_del_P64(kMers, key);
    }
}

/// Return the next k-mer in the k-mer set and update the index.
int64_t nextKMer(kh_P64_t *kMers, size_t &lastIndex) {
    for (size_t i = kh_begin(kMers) + lastIndex; i != kh_end(kMers); ++i, ++lastIndex) {
        if (!kh_exist(kMers, i)) continue;
        return kh_val(kMers, i);
    }
    // No more k-mers.
    return -1;
}



