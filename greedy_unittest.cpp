#include "greedy.cpp"

#include "gtest/gtest.h"

namespace {
    TEST(OverlapHamiltonianPathTest, OverlapHamiltonianPath) {
        struct TestCase {
            std::vector<KMer> kMers;
            std::vector<OverlapEdge> wantResult;
            int k;
        };
        std::vector<TestCase> tests = {
                {
                        std::vector<KMer>{KMer{"ACG"}, KMer{"TAC"}, KMer{"GGC"}},
                        std::vector<OverlapEdge>{OverlapEdge{1, 0, 2},OverlapEdge{0, 2, 1}},
                        3,
                },
        };

        for (auto t : tests) {
            std::vector<OverlapEdge> got = OverlapHamiltonianPath( t.kMers, t.k);
            EXPECT_EQ(t.wantResult.size(), got.size());
            for (int i = 0; i < t.wantResult.size(); ++i) {
                EXPECT_EQ(t.wantResult[i].firstIndex, got[i].firstIndex);
                EXPECT_EQ(t.wantResult[i].secondIndex, got[i].secondIndex);
                EXPECT_EQ(t.wantResult[i].overlapLength, got[i].overlapLength);
            }
        }
    }

    TEST(GreedyTest, Greedy) {
        std::vector<std::pair<KMerSet, std::vector<KMer>>> tests = {
                {KMerSet{"TACGT",  std::vector<bool> {1, 1, 1, 0, 0}, 3 }, {KMer{"CGT"}, KMer{"TAC"}, KMer{"ACG"}}},
                {KMerSet{"ACGTTT",  std::vector<bool> {1, 1, 0, 1, 0, 0}, 3 }, {KMer{"CGT"}, KMer{"TTT"}, KMer{"ACG"}}},
                {KMerSet{"TACTT",  std::vector<bool> {1, 1, 0, 0, 0}, 4 }, {KMer{"TACT"}, KMer{"ACTT"}}},
                {KMerSet{"TACTTAAGGAC",  std::vector<bool> {1, 1, 0, 0, 1, 0, 0, 1, 0, 0, 0}, 4 }, {KMer{"TACT"}, KMer{"ACTT"}, KMer{"GGAC"}, KMer{"TAAG"}}},
                {KMerSet{"GAAAAGTTTAAAGAC", std::vector<bool> {1,1, 0, 1,1,1,1,1,1,1,1,1,0,0,0}, 4}, {KMer{"AAGA"}, KMer{"TTAA"}, KMer{"TTTA"}, KMer{"AGAC"}, KMer{"GTTT"}, KMer{"AGTT"}, KMer{"AAGT"}, KMer{"TAAA"}, KMer{"AAAG"}, KMer{"AAAA"}, KMer{"GAAA"}}},
        };

        for (auto t : tests) {
            KMerSet got = Greedy(t.second);
            EXPECT_EQ(t.first.superstring, got.superstring);
            EXPECT_EQ(t.first.k, got.k);
            EXPECT_EQ(t.first.mask, got.mask);
        }
    }
}
