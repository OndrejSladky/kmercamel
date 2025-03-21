#pragma once
#include "../src/global.h"

#include <algorithm>

#include "kmer_types.h"

#include "gtest/gtest.h"
typedef unsigned char byte;
namespace {
    TEST(Global, PartialPreSort) {
        struct TestCase {
            std::vector<kmer_t> kMers;
            int k;
            std::vector<kmer_t> wantResult;
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

    TEST(Global, SuperstringFromPath) {
        struct TestCase {
            overlapPath path;
            std::vector<simplitig_t> kMers;
            int k;
            std::string wantResult;
            bool complements;
        };
        std::vector<TestCase> tests = {
                {
                        {{2, 0, (size_t)-1}, {1, 2, (byte)-1}},
                        std::vector<simplitig_t>{simplitig_from_string({"ACG"}), simplitig_from_string({"TAC"}), simplitig_from_string({"GGC"})},
                        3,
                        "TAcGgc",
                        false,
                },
                {
                        {{4, 3, 1, (size_t)-1, 5, (size_t)-1}, {1 ,2, 1, (byte)-1, 2, (byte)-1}},
                        std::vector<simplitig_t>{simplitig_from_string({"GCC"}), simplitig_from_string({"ACG"}), simplitig_from_string({"TAC"})},
                        3,
                        "GcCGta",
                        true,
                },
                {
                        {{3, 1, (size_t)-1, (size_t)-1}, {1, 1, (byte)-1, (byte)-1}},
                        std::vector<simplitig_t>{simplitig_from_string({"GCC"}), simplitig_from_string({"TACG"})},
                        3,
                        "GcCGta",
                        true,
                },
        };

        for (auto t : tests) {
            std::stringstream of;

            SuperstringFromPath(wrapper, kmer_t(0), t.path, t.kMers, of, nullptr, t.k, t.complements);

            EXPECT_EQ(t.wantResult, of.str());
        }
    }

    TEST(Global, OverlapHamiltonianPath) {
        struct TestCase {
            std::vector<simplitig_t> kMers;
            overlapPath wantResult;
            int k;
            bool complements;
            bool lower_bound;
        };
        std::vector<TestCase> tests = {
                {
                        {simplitig_from_string({"AT"})},
                        {{(size_t)-1, (size_t)-1}, {(byte)-1, (byte)-1}},
                        2,
                        true,
                        false,
                },
                {
                        {simplitig_from_string({"ACG"}), simplitig_from_string({"TAC"}), simplitig_from_string({"GGC"})},
                        {{2, 0, (size_t)-1}, {1, 2, (byte)-1}},
                        3,
                        false,
                        false,
                },
                {
                        {simplitig_from_string({"ACAA"}), simplitig_from_string({"ATTT"}), simplitig_from_string({"AACA"})},
                        {{4, 3, 0, 5, (size_t)-1, (size_t)-1}, {2, 2, 3, 3, (byte)-1, (byte)-1}},
                        4,
                        true,
                        false,
                },
                {
                        {simplitig_from_string({"ACG"}), simplitig_from_string({"CGT"}), simplitig_from_string({"TAA"})},
                        {{1, 2, 0}, {2, 1, 1}},
                        3,
                        false,
                        true,
                },
                {
                        {simplitig_from_string({"ACGT"}), simplitig_from_string({"TAA"})},
                        {{1, 0}, {1, 1}},
                        3,
                        false,
                        true,
                },
                {
                        {simplitig_from_string({"ACC"}), simplitig_from_string({"CGG"})},
                        {{3, 2, 1, 0}, {2, 2, 0, 0}},
                        3,
                        true,
                        true,
                },
        };

        for (auto t : tests) {
            overlapPath got = OverlapHamiltonianPath(wrapper, kmer_t(0), t.kMers, t.k, t.complements, t.lower_bound);
            EXPECT_EQ(t.wantResult.first, got.first);
            EXPECT_EQ(t.wantResult.second, got.second);
        }
    }

    TEST(Global, Global) {
        struct TestCase {
            std::string wantResult;
            int k;
            std::vector<simplitig_t> input;
            bool complements;
        };
        std::vector<TestCase> tests = {
                {"TACgt", 3, {simplitig_from_string({"CGT"}), simplitig_from_string({"TAC"}), simplitig_from_string({"ACG"})}, false},
                {"ACGT", 1, {simplitig_from_string({"ACGT"})}, false},
                {"ACgTtt", 3, {simplitig_from_string({"CGT"}), simplitig_from_string({"TTT"}), simplitig_from_string({"ACG"})}, false},
                {"ACgTtt", 3, {simplitig_from_string({"ACGT"}), simplitig_from_string({"TTT"})}, false},
                {"TActt", 4, {simplitig_from_string({"TACT"}), simplitig_from_string({"ACTT"})}, false},
                {"TActTaaGgac", 4, {simplitig_from_string({"TACT"}), simplitig_from_string({"ACTT"}), simplitig_from_string({"GGAC"}), simplitig_from_string({"TAAG"})}, false},
                {"TTtcttttttttttttttttttttttttttga", 31, {simplitig_from_string({"TTTCTTTTTTTTTTTTTTTTTTTTTTTTTTG"}), simplitig_from_string({"TTCTTTTTTTTTTTTTTTTTTTTTTTTTTGA"})}, false},
                {"AtTTgtt", 4, {simplitig_from_string({"ACAA"}), simplitig_from_string({"ATTT"}), simplitig_from_string({"AACA"})}, true},
        };

        for (auto &&t : tests) {
            std::stringstream of;

            Global(wrapper, kmer_t(0), t.input, of, nullptr, t.k, t.complements);

            EXPECT_EQ(t.wantResult, of.str());
        }
    }
}
