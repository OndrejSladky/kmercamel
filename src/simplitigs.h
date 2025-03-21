#pragma once

#include <vector>
#include "kmers.h"
#include "parser.h"
#include "local.h"
#include "khash_utils.h"

typedef std::vector<bool> simplitig_t;


simplitig_t simplitig_reverse_complement(simplitig_t forward) {
    simplitig_t reverse(forward.size());
    for (size_t i = 0; i < forward.size(); ++i) {
        size_t rc_position = (forward.size() - i - 1) ^ 1;
        reverse[rc_position] = 1 ^ forward[i];
    }
    return reverse;
}

char simplitig_at_index(simplitig_t &simplitig, size_t index) {
    return letters[simplitig[2 * index] * 2 + simplitig[2 * index + 1]];
}

inline size_t kmers_in_simplitig(simplitig_t &simplitig, int k) {
    return simplitig.size() / 2 - k + 1;
}

template <typename kmer_t>
kmer_t kmer_at_simplitig_index(kmer_t kmer_type, simplitig_t &simplitig, int k, size_t index) {
    kmer_t ret = 0;
    for (size_t i = 2 * index; i < 2 * (index + k); ++i) {
        ret <<= 1;
        ret |= simplitig[i];
    }
    return ret;
}


template <typename kmer_t>
kmer_t simplitig_first(kmer_t kmer_type, simplitig_t &simplitig, int k) {
    return kmer_at_simplitig_index(kmer_type, simplitig, k, 0);
}

template <typename kmer_t>
kmer_t simplitig_last(kmer_t kmer_type, simplitig_t &simplitig, int k) {
    return kmer_at_simplitig_index(kmer_type, simplitig, k, kmers_in_simplitig(simplitig, k) - 1); 
}

simplitig_t simplitig_from_string(std::string sequence) {
    simplitig_t ret(2 * sequence.size());
    for (size_t i = 0; i < sequence.size(); ++i) {
        uint8_t x = nucleotideToInt[(uint8_t)sequence[i]];
        assert(x <= 3);
        ret[2 * i] = x & 2;
        ret[2 * i + 1] = x & 1;
    }
    return ret;
}

std::vector<simplitig_t> simplitigs_from_fasta(std::string &path, int k) {
    std::vector<simplitig_t> simplitigs;
    
    gzFile fp = OpenFile(path);
    kseq_t *seq = kseq_init(fp);

    while (kseq_read(seq) >= 0) {
        simplitigs.push_back(simplitig_from_string(seq->seq.s));
    }

    kseq_destroy(seq);
    gzclose(fp);

    return simplitigs;
}

/// Find the next simplitig.
/// Also remove the used k-mers from kMers.
/// If complements are true, it is expected that kMers only contain one k-mer from a complementary pair.
template <typename kmer_t, typename kh_S_t, typename kh_wrapper_t>
simplitig_t next_simplitig(kh_S_t *kMers, kh_wrapper_t wrapper, kmer_t begin,  int k, bool complements) {
     // Maintain the first and last k-mer in the simplitig.
    kmer_t last = begin, first = begin;
    std::list<char> simplitig {};
    for (int i = 0; i < k - 1; ++i) {
        simplitig.emplace_back(NucleotideAtIndex(last, k, i));
    }
    eraseKMer(kMers, wrapper, last, k, complements);
    while (true) {
        auto extension = RightExtension(last, kMers, wrapper, k, 1, complements);
        kmer_t ext = extension.first;
        if (ext == kmer_t(-1)) {
            // No right extension found.
            break;
        } else {
            // Extend the simplitig to the right.
            eraseKMer(kMers, wrapper, extension.second, k, complements);
            simplitig.emplace_back(letters[(uint8_t)ext]);
            last = extension.second;
        }
    } 
    while(true) {
        auto extension = LeftExtension(first, kMers,  wrapper, k, 1, complements);
        kmer_t ext = extension.first;
        if (ext == kmer_t(-1)) {
            // No left extension found.
            break;
        } else {
            // Extend the simplitig to the left.
            eraseKMer(kMers, wrapper, extension.second, k, complements);
            simplitig.emplace_front(letters[(uint8_t)ext]);
            first = extension.second;
        }
    }

    return simplitig_from_string(std::string(simplitig.begin(), simplitig.end()));
}

template <typename kmer_t, typename kh_S_t, typename kh_wrapper_t>
std::vector<simplitig_t> get_simplitigs(kh_S_t *kMers, kh_wrapper_t wrapper, kmer_t _, int k, bool complements) {
    size_t lastIndex = 0;
    std::vector<simplitig_t> simplitigs;
    while(true) {
        kmer_t begin = nextKMer(kMers, _, lastIndex);
        // No more k-mers.
        if (begin == kmer_t(-1)) return simplitigs;
        simplitigs.emplace_back(next_simplitig(kMers, wrapper, begin, k, complements));
    }
    return simplitigs;
}
