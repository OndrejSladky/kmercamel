#ifndef KMERS_CPP
#define KMERS_CPP
#include <string>
#include <iostream>

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

long long KMerToNumber(KMer &kMer) {
    long long ret = 0;
    for (char c : kMer.value) {
        ret <<= 2;
        ret |= NucleotideToInt(c);
    }
    return ret;
}

long long BitPrefix(long long kMer, int k, int d) {
    long long mask = -1LL ^ ((1LL << ((k - d) << 1LL)) - 1LL);
    return (kMer & mask) >> ((k - d) << 1LL);
}

long long BitSuffix(long long kMer, int d) {
    return kMer & ((1LL << (d << 1LL)) - 1LL);
}

long long ReverseComplement(long long kMer, int k) {
    long long ans = 0;
    for (int i = 0; i < k; ++i) {
        ans <<= 2;
        ans |= 3 ^ (kMer & 3);
        kMer >>= 2;
    }
    return ans;
}

const char letters[4] {'A', 'C', 'G', 'T'};

/// Convert the encoded KMer representation to string.
std::string NumberToKMer(long long encoded, int length) {
    std::string ret = "";
    for (int i = 0; i < length; ++i) {
        // The last two bits correspond to one nucleotide.
        ret = letters[encoded & 3] + ret;
        // Move to the next letter.
        encoded >>= 2;
    }
    return ret;
}


#endif //KMERS_CPP