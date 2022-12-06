#pragma once
#include "greedy.h"

#include "gtest/gtest.h"

namespace {
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
                {
                        {KMerToNumber({"ACAA"}), KMerToNumber({"ATTT"}),KMerToNumber({"CCCC"}), KMerToNumber({"AACA"}), KMerToNumber({"TTGT"}), KMerToNumber({"AAAT"}),KMerToNumber({"GGGG"}), KMerToNumber({"TGTT"})},
                        std::vector<OverlapEdge>{{3, 0, 3},{4, 7, 3},{0, 5, 2},{1, 4, 2}, {2, 1, 0},{5, 6, 0}},
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
            KMerSet wantResult;
            std::vector<int64_t> input;
            bool complements;
        };
        std::vector<TestCase> tests = {
                {KMerSet{"TACGT", std::vector<bool> {1, 1, 1, 0, 0}, 3 }, {KMerToNumber({"CGT"}), KMerToNumber({"TAC"}), KMerToNumber({"ACG"})}, false},
                {KMerSet{"TAGC", std::vector<bool> {1, 0, 1, 0}, 2 }, {KMerToNumber({"TA"}), KMerToNumber({"GC"}), }, false},
                {KMerSet{"ACGTTT", std::vector<bool> {1, 1, 0, 1, 0, 0}, 3 }, {KMerToNumber({"CGT"}), KMerToNumber({"TTT"}), KMerToNumber({"ACG"})}, false},
                {KMerSet{"TACTT", std::vector<bool> {1, 1, 0, 0, 0}, 4 }, {KMerToNumber({"TACT"}), KMerToNumber({"ACTT"})}, false},
                {KMerSet{"TACTTAAGGAC",  std::vector<bool> {1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0}, 4 }, {KMerToNumber({"TACT"}), KMerToNumber({"ACTT"}), KMerToNumber({"GGAC"}), KMerToNumber({"TAAG"})}, false},
                {KMerSet{"GAAAAGTTTAAAGAC", std::vector<bool> {1,1, 0, 1,1,1,1,1,1,1,1,1,0,0,0}, 4}, {KMerToNumber({"AAGA"}), KMerToNumber({"TTAA"}), KMerToNumber({"TTTA"}), KMerToNumber({"AGAC"}), KMerToNumber({"GTTT"}), KMerToNumber({"AGTT"}), KMerToNumber({"AAGT"}), KMerToNumber({"TAAA"}), KMerToNumber({"AAAG"}), KMerToNumber({"AAAA"}), KMerToNumber({"GAAA"})}, false},
                {KMerSet{"TTTCTTTTTTTTTTTTTTTTTTTTTTTTTTGA", std::vector<bool> {1,1,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0}, 31}, {KMerToNumber({"TTTCTTTTTTTTTTTTTTTTTTTTTTTTTTG"}), KMerToNumber({"TTCTTTTTTTTTTTTTTTTTTTTTTTTTTGA"})}, false},
                {KMerSet{"CCCCATTTGTT", {1,0,0,0,1,0,1,1,0,0,0,}, 4}, {KMerToNumber({"ACAA"}), KMerToNumber({"ATTT"}), KMerToNumber({"CCCC"}), KMerToNumber({"AACA"})}, true},
        };

        for (auto &&t : tests) {
            KMerSet got = Greedy(t.input, t.wantResult.k, t.complements);
            EXPECT_EQ(t.wantResult.superstring, got.superstring);
            EXPECT_EQ(t.wantResult.k, got.k);
            EXPECT_EQ(t.wantResult.mask, got.mask);
        }
    }
}
