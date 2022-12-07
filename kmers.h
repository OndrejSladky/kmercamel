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
int64_t KMerToNumber(const KMer &kMer) {
    int64_t ret = 0;
    for (char c : kMer.value) {
        ret <<= 2;
        ret |= NucleotideToInt(c);
    }
    return ret;
}

/// Compute the prefix of size d of the given k-mer.
int64_t BitPrefix(int64_t kMer, int k, int d) {
    return kMer >> ((k - d) << 1LL);
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

/// Return the complementary nucleotide for the given one.
char ComplementaryNucleotide(const char nucleotide) {
    if (nucleotide == 'A') return 'T';
    else if (nucleotide == 'T') return 'A';
    else if (nucleotide == 'G') return 'C';
    else if (nucleotide == 'C') return 'G';
    throw std::invalid_argument("cannot find complementary nucleotide for letter " + std::string(1, nucleotide));
}

/// Compute the reverse complement of the given k-mer.
KMer ReverseComplement(const KMer &kMer) {
    KMer ans;
    ans.value.reserve(kMer.length());
    for (int i = (int)kMer.length() - 1; i >= 0; --i) {
        ans.value += ComplementaryNucleotide(kMer.value[i]);
    }
    return ans;
}

const char letters[4] {'A', 'C', 'G', 'T'};

/// Convert the encoded KMer representation to string.
std::string NumberToKMer(int64_t encoded, int length) {
    std::string ret;
    ret.reserve(length);
    for (int i = 0; i < length; ++i) {
        // The last two bits correspond to one nucleotide.
        ret[length - i -1] = letters[encoded & 3];
        // Move to the next letter.
        encoded >>= 2;
    }
    return ret;
}
