#pragma once

#include <vector>
#include <cassert>
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

template <bool complements, typename kmer_t, typename kh_S_t, typename kh_wrapper_t>
inline kmer_t simplitig_right_rev_extension(kmer_t &forward, kmer_t &backward, kh_S_t *kMers, kh_wrapper_t &wrapper, int k) {
    bool forward_direction = true;
    if constexpr (!complements) {
        forward_direction = (backward == kmer_t(-1));
    }
    kmer_t mask = (kmer_t(1) << (k << 1)) - kmer_t(1);
    kmer_t forward_base = (forward << 2) & mask;
    kmer_t backward_base = backward >> 2;
    for (kmer_t ext = 0; ext < kmer_t(4); ++ext) {
        kmer_t next_forward = forward_base | (kmer_t(3) ^ ext);
        kmer_t next_backward = backward_base | (ext << (((k - 1) << 1)));
        kmer_t canonical;
        if constexpr (complements) {
            canonical = (next_forward <= next_backward) ? next_forward : next_backward;
        } else {
            canonical = (forward_direction) ? next_forward : next_backward;
        }
        auto key = wrapper.kh_get_from_set(kMers, canonical);
        if(key != kh_end(kMers)) {
            if constexpr (complements) {
                backward = next_backward;
                forward = next_forward;
            } else {
                if (forward_direction) forward = next_forward;
                else backward = next_backward;
            }
            wrapper.kh_del_from_set(kMers, key);
            return ext;
        }
    }
    return -1;
}

/// Find the next simplitig.
/// Also remove the used k-mers from kMers.
/// If complements are true, it is expected that kMers only contain one k-mer from a complementary pair.
template <bool complements, typename kmer_t, typename kh_S_t, typename kh_wrapper_t>
simplitig_t next_simplitig(kh_S_t *kMers, kh_wrapper_t wrapper, kmer_t begin, int k) {
     // Maintain the first and last k-mer in the simplitig.
    kmer_t last = begin, first = begin;
    kmer_t last_complement = ReverseComplement(begin, k);
    kmer_t first_complement = last_complement;
    if constexpr (!complements) {
        first_complement = last_complement = kmer_t(-1);
    }
    std::list<char> simplitig {};
    for (int i = 0; i < k; ++i) {
        simplitig.emplace_back(NucleotideAtIndex(last, k, i));
    }
    eraseKMer(kMers, wrapper, last, k, complements);
    while (true) {
        kmer_t ext = simplitig_right_rev_extension<complements>(last, last_complement, kMers, wrapper, k);
        if (ext == kmer_t(-1)) {
            // No right extension found.
            break;
        } else {
            // Extend the simplitig to the right.
            simplitig.emplace_back(letters[3 ^ (uint8_t)ext]);
        }
    } 
    while(true) {
        kmer_t ext = simplitig_right_rev_extension<complements>(first_complement, first, kMers,  wrapper, k);
        if (ext == kmer_t(-1)) {
            // No left extension found.
            break;
        } else {
            // Extend the simplitig to the left.
            simplitig.emplace_front(letters[(uint8_t)ext]);
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
        simplitig_t next;
        if (complements) next = next_simplitig<true>(kMers, wrapper, begin, k);
        else next = next_simplitig<false>(kMers, wrapper, begin, k);
        simplitigs.emplace_back(next);
    }
    return simplitigs;
}
