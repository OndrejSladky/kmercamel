#include "greedy.cpp"

#include "gtest/gtest.h"

namespace {
    TEST(OverlapHamiltonianPathTest, OverlapHamiltonianPath) {
        struct TestCase {
            std::vector<KMer> kMers;
            std::vector<OverlapEdge> wantResult;
            int k;
            bool complements;
        };
        std::vector<TestCase> tests = {
                {
                        {KMer{"AT"},  },
                        std::vector<OverlapEdge>{},
                        2,
                        true,
                },
                {
                        std::vector<KMer>{KMer{"ACG"}, KMer{"TAC"}, KMer{"GGC"}},
                        std::vector<OverlapEdge>{OverlapEdge{1, 0, 2},OverlapEdge{0, 2, 1}},
                        3,
                        false,
                },
                {
                        {KMer{"ACAA"}, KMer{"ATTT"}, KMer{"AACA"}},
                        std::vector<OverlapEdge>{OverlapEdge{2, 0, 3},OverlapEdge{0, 4, 2}},
                        4,
                        true,
                },
                {
                        {KMer{"ACAA"}, KMer{"ATTT"}, KMer{"CCCC"}, KMer{"AACA"}},
                        std::vector<OverlapEdge>{OverlapEdge{3, 0, 3},OverlapEdge{0, 5, 2}, {2, 3, 0}},
                        4,
                        true,
                },
        };

        for (auto t : tests) {
            std::vector<OverlapEdge> got = OverlapHamiltonianPath( t.kMers, t.k, t.complements);
            EXPECT_EQ(t.wantResult.size(), got.size());
            for (int i = 0; i < t.wantResult.size(); ++i) {
                EXPECT_EQ(t.wantResult[i].firstIndex, got[i].firstIndex);
                EXPECT_EQ(t.wantResult[i].secondIndex, got[i].secondIndex);
                EXPECT_EQ(t.wantResult[i].overlapLength, got[i].overlapLength);
            }
        }
    }

    TEST(GreedyTest, Greedy) {
        struct TestCase {
            KMerSet wantResult;
            std::vector<KMer> input;
            bool complements;
        };
        std::vector<TestCase> tests = {
                {KMerSet{"TACGT", std::vector<bool> {1, 1, 1, 0, 0}, 3 }, {KMer{"CGT"}, KMer{"TAC"}, KMer{"ACG"}}, false},
                {KMerSet{"TAGC", std::vector<bool> {1, 0, 1, 0}, 2 }, {KMer{"TA"}, KMer{"GC"}, }, false},
                {KMerSet{"ACGTTT", std::vector<bool> {1, 1, 0, 1, 0, 0}, 3 }, {KMer{"CGT"}, KMer{"TTT"}, KMer{"ACG"}}, false},
                {KMerSet{"TACTT", std::vector<bool> {1, 1, 0, 0, 0}, 4 }, {KMer{"TACT"}, KMer{"ACTT"}}, false},
                {KMerSet{"TACTTAAGGAC",  std::vector<bool> {1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0}, 4 }, {KMer{"TACT"}, KMer{"ACTT"}, KMer{"GGAC"}, KMer{"TAAG"}}, false},
                {KMerSet{"GAAAAGTTTAAAGAC", std::vector<bool> {1,1, 0, 1,1,1,1,1,1,1,1,1,0,0,0}, 4}, {KMer{"AAGA"}, KMer{"TTAA"}, KMer{"TTTA"}, KMer{"AGAC"}, KMer{"GTTT"}, KMer{"AGTT"}, KMer{"AAGT"}, KMer{"TAAA"}, KMer{"AAAG"}, KMer{"AAAA"}, KMer{"GAAA"}}, false},
                {KMerSet{"TTTCTTTTTTTTTTTTTTTTTTTTTTTTTTGA", std::vector<bool> {1,1,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0,  0,0,0,0}, 31}, {KMer{"TTTCTTTTTTTTTTTTTTTTTTTTTTTTTTG"}, KMer{"TTCTTTTTTTTTTTTTTTTTTTTTTTTTTGA"}}, false},
                {KMerSet{"CCCCAACAAAT", {1,0,0,0,1,1,0,1,0,0,0,}, 4}, {KMer{"ACAA"}, KMer{"ATTT"}, KMer{"CCCC"}, KMer{"AACA"}}, true},
        };

        for (auto t : tests) {
            KMerSet got = Greedy(t.input, t.complements);
            EXPECT_EQ(t.wantResult.superstring, got.superstring);
            EXPECT_EQ(t.wantResult.k, got.k);
            EXPECT_EQ(t.wantResult.mask, got.mask);
        }
    }
}
