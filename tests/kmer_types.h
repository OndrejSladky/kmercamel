#pragma once

#include <cstdint>
#include "../src/khash_utils.h"


#ifdef LARGE_KMERS
    typedef __uint128_t kmer_t;
    typedef kmer_dict128_t kh_wrapper;
    typedef kh_S128_t kh_S_t;
    typedef kh_P128_t kh_P_t;
#else // LARGE_KMERS
    typedef uint64_t kmer_t;
    typedef kmer_dict64_t kh_wrapper;
    typedef kh_S64_t kh_S_t;
    typedef kh_P64_t kh_P_t;
#endif // LARGE_KMERS

// To be passed to functions.
kh_wrapper wrapper;