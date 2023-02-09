#pragma once
#include "greedy.h"

#include "gtest/gtest.h"

namespace {
    TEST(SuperstringFromPathTest, Number) {
        struct TestCase {
            std::vector<OverlapEdge> path;
            std::vector<int64_t> kMers;
            int k;
            std::string wantResult;
        };
        std::vector<TestCase> tests = {
                {
                        std::vector<OverlapEdge>{OverlapEdge{1, 0, 2},OverlapEdge{0, 2, 1}},
                        std::vector<int64_t>{KMerToNumber({"ACG"}), KMerToNumber({"TAC"}), KMerToNumber({"GGC"})},
                        3,
                        "TAcGgc",
                },
                {
                        std::vector<OverlapEdge>{OverlapEdge{2, 1, 2},OverlapEdge{1, 3, 1}},
                        std::vector<int64_t>{KMerToNumber({"GCC"}), KMerToNumber({"ACG"}), KMerToNumber({"TAC"}),
                                          KMerToNumber({"GGC"}), KMerToNumber({"CGT"}), KMerToNumber({"GTA"})},
                        3,
                        "TAcGgc",
                },
        };

        for (auto t : tests) {
            std::stringstream of;

            SuperstringFromPath(t.path, t.kMers, of, t.k);

            EXPECT_EQ(t.wantResult, of.str());
        }
    }

    TEST(OverlapHamiltonianPathTest, OverlapHamiltonianPath) {
        struct TestCase {
            std::vector<int64_t> kMers;
            std::vector<OverlapEdge> wantResult;
            int k;
            bool complements;
        };
        std::vector<TestCase> tests = {
                {
                        {KMerToNumber({"AT"}), KMerToNumber({"AT"}) },
                        std::vector<OverlapEdge>{},
                        2,
                        true,
                },
                {
                        {KMerToNumber({"ACG"}), KMerToNumber({"TAC"}), KMerToNumber({"GGC"})},
                        std::vector<OverlapEdge>{OverlapEdge{1, 0, 2},OverlapEdge{0, 2, 1}},
                        3,
                        false,
                },
                {
                        {KMerToNumber({"ACAA"}), KMerToNumber({"ATTT"}), KMerToNumber({"AACA"}), KMerToNumber({"TTGT"}), KMerToNumber({"AAAT"}), KMerToNumber({"TGTT"})},
                        std::vector<OverlapEdge>{{2, 0, 3},{3, 5, 3},{0, 4, 2},{1, 3, 2}},
                        4,
                        true,
                },
        };

        for (auto t : tests) {
            std::vector<OverlapEdge> got = OverlapHamiltonianPath( t.kMers, t.k, t.complements);
            EXPECT_EQ(t.wantResult.size(), got.size());
            for (size_t i = 0; i < t.wantResult.size(); ++i) {
                EXPECT_EQ(t.wantResult[i].firstIndex, got[i].firstIndex);
                EXPECT_EQ(t.wantResult[i].secondIndex, got[i].secondIndex);
                EXPECT_EQ(t.wantResult[i].overlapLength, got[i].overlapLength);
            }
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
                {"TaGc", 2, {KMerToNumber({"TA"}), KMerToNumber({"GC"}), }, false},
                {"ACgTtt", 3, {KMerToNumber({"CGT"}), KMerToNumber({"TTT"}), KMerToNumber({"ACG"})}, false},
                {"TActt", 4, {KMerToNumber({"TACT"}), KMerToNumber({"ACTT"})}, false},
                {"TActTaaGgac", 4, {KMerToNumber({"TACT"}), KMerToNumber({"ACTT"}), KMerToNumber({"GGAC"}), KMerToNumber({"TAAG"})}, false},
                {"TTtcttttttttttttttttttttttttttga", 31, {KMerToNumber({"TTTCTTTTTTTTTTTTTTTTTTTTTTTTTTG"}), KMerToNumber({"TTCTTTTTTTTTTTTTTTTTTTTTTTTTTGA"})}, false},
                {"AtTTgttGggg", 4, {KMerToNumber({"ACAA"}), KMerToNumber({"ATTT"}), KMerToNumber({"CCCC"}), KMerToNumber({"AACA"})}, true},
        };

        for (auto &&t : tests) {
            std::stringstream of;

            Greedy(t.input, of, t.k, t.complements);

            EXPECT_EQ(t.wantResult, of.str());
        }
    }
}
