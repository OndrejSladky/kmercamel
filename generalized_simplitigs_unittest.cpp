#include "generalized_simplitigs.cpp"

#include "gtest/gtest.h"

namespace {
    TEST(NumberToKMerTest, NumberToKMer) {
        struct TestCase {
            long long encoded;
            int d;
            std::string wantResult;
        };
        std::vector<TestCase> tests = {
                {0b1001LL, 2, "GC"},
                {0b1011LL, 3, "AGT"},
                {0b111LL, 1, "T"},
        };

        for (auto t: tests) {
            std::string gotResult = NumberToKMer(t.encoded, t.d);

            EXPECT_EQ(t.wantResult, gotResult);
        }
    }

    TEST(ExtensionTest, Extension) {
        struct TestCase {
            std::string simplitig;
            std::unordered_set<std::string> kMers;
            int k;
            int d;
            bool fromRight;
            std::string wantExt;
            std::string wantNext;
        };
        std::vector<TestCase> tests = {
                {"GACT", {"TCC", "CTA", "ACT", "CCT"}, 3, 1, true,   "A", "CTA"},
                {"TTT", {"TCC",  "ACT", "CCT"}, 3, 2, true,   "CC", "TCC"},
                {"TACT", {"TCC",  "ACT", "CCT"}, 3, 1, true,   "", ""},
                {"TACG", {"TCC", "CTA", "ACT", "CCT"}, 3, 1, false,   "C", "CTA"},
        };

        for (auto t: tests) {
            auto got = Extension(t.simplitig, t.kMers, t.k, t.d, t.fromRight);
            auto gotExt = got.first;
            auto gotNext = got.second;
            EXPECT_EQ(t.wantNext, gotNext);
            EXPECT_EQ(t.wantExt, gotExt);
        }
    }

    TEST(NextGeneralizedSimplitigTest, NextGeneralizedSimplitig) {
        struct TestCase {
            std::string superstring;
            std::unordered_set<std::string> kMers;
            std::vector<bool> mask;
            int k;
            int d_max;
            std::string wantSuperstring;
            std::unordered_set<std::string> wantKMers;
            std::vector<bool> wantMask;
        };
        std::vector<TestCase> tests = {
                // Behavior of the unordered_set.begin() is not specified. It may happen that it takes "ATTT" first and then this test fails.
                {"GGGG", {"ACAA", "ATTT", "AACA"}, {1,0,0, 0}, 4, 2,
                 "GGGGAACAA", {"ATTT"}, {1,0,0,0,1,1,0,0, 0}},
        };

        for (auto t: tests) {
            NextGeneralizedSimplitig(t.kMers, t.superstring, t.mask, t.k, t.d_max);

            EXPECT_EQ(t.wantSuperstring, t.superstring);
            EXPECT_EQ(t.wantMask, t.mask);
            EXPECT_EQ(t.wantKMers, t.kMers);
        }
    }

    TEST(GreedyGeneralizedSimplitigsTest, GreedyGeneralizedSimplitigs) {
        struct TestCase {
            std::vector<KMer> kMers;
            int k;
            int d_max;
            std::string wantSuperstring;
            std::vector<bool> wantMask;
        };
        std::vector<TestCase> tests = {
                // Behavior of the unordered_set.begin() is not specified. It may happen that it takes "CCCC" first and then this test fails.
                { {KMer{"ACAA"}, KMer{"ATTT"}, KMer{"CCCC"}, KMer{"AACA"}}, 4, 3,
                 "AACAATTTCCCC", {1,1,0,0,1, 0,0,0, 1, 0,0,0}},
        };

        for (auto t: tests) {
            auto gotResult = GreedyGeneralizedSimplitigs(t.kMers, t.k, t.d_max);

            EXPECT_EQ(t.wantSuperstring, gotResult.superstring);
            EXPECT_EQ(t.wantMask, gotResult.mask);
            EXPECT_EQ(t.k, gotResult.k);
        }
    }


}
