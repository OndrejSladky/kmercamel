#pragma once
#include "greedy.h"

#include <algorithm>

#include "gtest/gtest.h"
typedef unsigned char byte;
namespace {
    TEST(PartialPreSort, PartialPreSort) {
        struct TestCase {
            std::vector<int64_t> kMers;
            int k;
            std::vector<int64_t> wantResult;
        };
        std::vector<TestCase> tests = {
                {
                        {KMerToNumber({"GTA"}), KMerToNumber({"TAC"}), KMerToNumber({"GGC"})},
                        3,
                        {KMerToNumber({"GGC"}), KMerToNumber({"GTA"}), KMerToNumber({"TAC"})},
                },
                {
                        {KMerToNumber({"TTTTTTTTTTTTT"}),KMerToNumber({"AAAAAAAAAAAAA"}),  KMerToNumber({"GCGCGCGCGCGCG"})},
                        13,
                        {KMerToNumber({"AAAAAAAAAAAAA"}), KMerToNumber({"GCGCGCGCGCGCG"}), KMerToNumber({"TTTTTTTTTTTTT"})},
                },
                {
                    // Sorting is done only based on the first 11 chars.
                        {KMerToNumber({"AAAAAAAAAAAAT"}),KMerToNumber({"AAAAAAAAAAAAA"})},
                        13,
                        {KMerToNumber({"AAAAAAAAAAAAT"}),KMerToNumber({"AAAAAAAAAAAAA"})},
                },
        };

        for (auto t : tests) {
            PartialPreSort(t.kMers, t.k);

            EXPECT_EQ(t.wantResult, t.kMers);
        }
    }

    TEST(SuperstringFromPathTest, Number) {
        struct TestCase {
            overlapPath path;
            std::vector<int64_t> kMers;
            int k;
            std::string wantResult;
            bool complements;
        };
        std::vector<TestCase> tests = {
                {
                        {{2, 0, (size_t)-1}, {1, 2, (byte)-1}},
                        std::vector<int64_t>{KMerToNumber({"ACG"}), KMerToNumber({"TAC"}), KMerToNumber({"GGC"})},
                        3,
                        "TAcGgc",
                        false,
                },
                {
                        {{4, 3, 1, (size_t)-1, 5, (size_t)-1}, {1 ,2, 1, (byte)-1, 2, (byte)-1}},
                        std::vector<int64_t>{KMerToNumber({"GCC"}), KMerToNumber({"ACG"}), KMerToNumber({"TAC"})},
                        3,
                        "GcCGta",
                        true,
                },
        };

        for (auto t : tests) {
            std::stringstream of;

            SuperstringFromPath(t.path, t.kMers, of, t.k, t.complements);

            EXPECT_EQ(t.wantResult, of.str());
        }
    }

    TEST(OverlapHamiltonianPathTest, OverlapHamiltonianPath) {
        struct TestCase {
            std::vector<int64_t> kMers;
            overlapPath wantResult;
            int k;
            bool complements;
        };
        std::vector<TestCase> tests = {
                {
                        {KMerToNumber({"AT"})},
                        {{(size_t)-1, (size_t)-1}, {(byte)-1, (byte)-1}},
                        2,
                        true,
                },
                {
                        {KMerToNumber({"ACG"}), KMerToNumber({"TAC"}), KMerToNumber({"GGC"})},
                        {{2, 0, (size_t)-1}, {1, 2, (byte)-1}},
                        3,
                        false,
                },
                {
                        {KMerToNumber({"ACAA"}), KMerToNumber({"ATTT"}), KMerToNumber({"AACA"})},
                        {{4, 3, 0, 5, (size_t)-1, (size_t)-1}, {2, 2, 3, 3, (byte)-1, (byte)-1}},
                        4,
                        true,
                },
        };

        for (auto t : tests) {
            overlapPath got = OverlapHamiltonianPath( t.kMers, t.k, t.complements);
            EXPECT_EQ(t.wantResult.first, got.first);
            EXPECT_EQ(t.wantResult.second, got.second);
        }
    }

    TEST(GreedyTest, Greedy) {
        struct TestCase {
            std::string wantResult;
            int k;
            std::vector<int64_t> input;
            bool complements;
        };
        std::vector<TestCase> tests = {
                {"TACgt", 3, {KMerToNumber({"CGT"}), KMerToNumber({"TAC"}), KMerToNumber({"ACG"})}, false},
                {"ACgTtt", 3, {KMerToNumber({"CGT"}), KMerToNumber({"TTT"}), KMerToNumber({"ACG"})}, false},
                {"TActt", 4, {KMerToNumber({"TACT"}), KMerToNumber({"ACTT"})}, false},
                {"TActTaaGgac", 4, {KMerToNumber({"TACT"}), KMerToNumber({"ACTT"}), KMerToNumber({"GGAC"}), KMerToNumber({"TAAG"})}, false},
                {"TTtcttttttttttttttttttttttttttga", 31, {KMerToNumber({"TTTCTTTTTTTTTTTTTTTTTTTTTTTTTTG"}), KMerToNumber({"TTCTTTTTTTTTTTTTTTTTTTTTTTTTTGA"})}, false},
                {"AtTTgtt", 4, {KMerToNumber({"ACAA"}), KMerToNumber({"ATTT"}), KMerToNumber({"AACA"})}, true},
        };

        for (auto &&t : tests) {
            std::stringstream of;

            Greedy(t.input, of, t.k, t.complements);

            EXPECT_EQ(t.wantResult, of.str());
        }
    }
}
