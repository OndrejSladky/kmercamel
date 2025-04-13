#pragma once

#include <vector>
#include <cassert>
#include <cmath>
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

template <typename kmer_t, typename kh_S_t, typename kh_wrapper_t>
void fill_kmers(kh_S_t kmers, kh_wrapper_t wrapper, kmer_t kmer_type, std::vector<simplitig_t> simplitigs, int k) {
    for (auto simplitig : simplitigs) {
        auto simplitig_kmers = kmers_in_simplitig(simplitig, k);
        for (size_t i = 0; i < simplitig_kmers; ++i) {
            auto kmer = kmer_at_simplitig_index(kmer_type, simplitig, k, i);
            int ret;
            wrapper.kh_put_to_set(kmers, kmer, &ret);
        }
    }
}

template <typename kmer_t>
std::vector<kmer_t> simplitigs_to_kmer_vec(kmer_t kmer_type, std::vector<simplitig_t> simplitigs, int k, size_t length_estimate) {
    std::vector<kmer_t> result;
    result.reserve(length_estimate);
    for (auto simplitig : simplitigs) {
        auto simplitig_kmers = kmers_in_simplitig(simplitig, k);
        for (size_t i = 0; i < simplitig_kmers; ++i) {
            auto kmer = kmer_at_simplitig_index(kmer_type, simplitig, k, i);
            result.push_back(kmer);
        }
    }
    return result;
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

std::vector<simplitig_t> simplitigs_from_fasta(std::string &path) {
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
inline uint8_t simplitig_right_rev_extension(kmer_t &forward, kmer_t &backward, kh_S_t *kMers, kh_wrapper_t &wrapper, int k) {
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
    // Assumes that the largest simplitig is in order of sqrt(n) to save some time resizing the vector.
    size_t size_estimate = std::sqrt(kh_size(kMers)) * 2;
    simplitig_t simplitig_front;
    simplitig_t simplitig_back(2 * k);
    simplitig_front.reserve(size_estimate);
    simplitig_back.reserve(2 * k + size_estimate);
    for (int i = 0; i < 2 * k; ++i) {
        simplitig_back[2 * k - i - 1] = last & (kmer_t(1) << i);
    }
    eraseKMer(kMers, wrapper, last, k, complements);
    while (true) {
        uint8_t ext = simplitig_right_rev_extension<complements>(last, last_complement, kMers, wrapper, k);
        if (ext == uint8_t(-1)) {
            // No right extension found.
            break;
        } else {
            // Extend the simplitig to the right.

            simplitig_back.push_back(!(ext & 2));
            simplitig_back.push_back(!(ext & 1));
        }
    } 
    while(true) {
        uint8_t ext = simplitig_right_rev_extension<complements>(first_complement, first, kMers,  wrapper, k);
        if (ext == uint8_t(-1)) {
            // No left extension found.
            break;
        } else {
            // Extend the simplitig to the left.
            simplitig_front.push_back(ext & 1);
            simplitig_front.push_back(ext & 2);
        }
    }

    std::reverse(simplitig_front.begin(), simplitig_front.end());
    simplitig_front.reserve(simplitig_front.size() + simplitig_back.size());
    for (bool b : simplitig_back) simplitig_front.push_back(b);
    return simplitig_front;
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
