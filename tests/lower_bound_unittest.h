#pragma once

#include "kmer_types.h"

#include "../src/lower_bound.h"

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
}