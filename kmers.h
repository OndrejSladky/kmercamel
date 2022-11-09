#pragma once
#include <string>
#include <iostream>
#include <cstdint>

#include "models.h"

/// Convert the given basic nucleotide to int so it can be used for indexing in AC.
int NucleotideToInt (char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: throw std::invalid_argument("cannot convert letter " + std::string(1, c) + "to int.");
    }
}

/// Convert the given k-mer to its representation as a number.
int64_t KMerToNumber(KMer &kMer) {
    int64_t ret = 0;
    for (char c : kMer.value) {
        ret <<= 2;
        ret |= NucleotideToInt(c);
    }
    return ret;
}

/// Compute the prefix of size d of the given k-mer.
int64_t BitPrefix(int64_t kMer, int k, int d) {
    int64_t mask = -1LL ^ ((1LL << ((k - d) << 1LL)) - 1LL);
    return (kMer & mask) >> ((k - d) << 1LL);
}

/// Compute the suffix of size d of the given k-mer.
int64_t BitSuffix(int64_t kMer, int d) {
    return kMer & ((1LL << (d << 1LL)) - 1LL);
}

/// Compute the reverse complement of the given k-mer.
int64_t ReverseComplement(int64_t kMer, int k) {
    int64_t ans = 0;
    for (int i = 0; i < k; ++i) {
        ans <<= 2;
        ans |= 3 ^ (kMer & 3);
        kMer >>= 2;
    }
    return ans;
}

const char letters[4] {'A', 'C', 'G', 'T'};

/// Convert the encoded KMer representation to string.
std::string NumberToKMer(int64_t encoded, int length) {
    std::string ret = "";
    for (int i = 0; i < length; ++i) {
        // The last two bits correspond to one nucleotide.
        ret = letters[encoded & 3] + ret;
        // Move to the next letter.
        encoded >>= 2;
    }
    return ret;
}
