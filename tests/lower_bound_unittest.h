#pragma once

#include "kmer_types.h"

#include "../src/lower_bound.h"
#include "../src/simplitigs.h"

#include "gtest/gtest.h"

namespace {
    TEST(LowerBound, LowerBoundLength) {
        struct TestCase {
            std::vector<simplitig_t> kMers;
            int k;
            bool complements;
            size_t wantResult;
        };
        std::vector<TestCase> tests = {
                {
                        {simplitig_from_string({"ACG"}), simplitig_from_string({"CGT"}), simplitig_from_string({"TAA"})},
                        3, false, 5
                },
                {
                        {simplitig_from_string({"ACGT"}), simplitig_from_string({"TAA"})},
                        3, false, 5
                },
                {
                        {simplitig_from_string({"AC"}), simplitig_from_string({"GG"})},
                        2, true, 3
                },
                {
                        {simplitig_from_string({"AAACCC"}), simplitig_from_string({"CCCAAA"}), simplitig_from_string({"TGGGGT"})},
                        6, false, 11
                },
        };

        for (auto &t : tests) {
            auto gotResult = LowerBoundLength(wrapper, kmer_t(0), t.kMers, t.k, t.complements);

            ASSERT_EQ(t.wantResult, gotResult);
        }
    }

    TEST(LowerBound, LowerBoundLengthSparseMatchesDenseForSingletonSimplitigs) {
        std::vector<simplitig_t> simplitigs = {
                simplitig_from_string({"ACG"}), simplitig_from_string({"CGT"}), simplitig_from_string({"TAA"})};
        constexpr int k = 3;
        auto kMerVec = simplitigs_to_kmer_vec(kmer_t(0), simplitigs, k, 16);
        size_t sparse = LowerBoundLengthSparse(wrapper, kMerVec, k, false);
        size_t dense = LowerBoundLength(wrapper, kmer_t(0), simplitigs, k, false);
        ASSERT_EQ(dense, sparse);
    }
}