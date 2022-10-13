#include "greedy.cpp"

#include "gtest/gtest.h"

namespace {
    TEST(SuffixTest, Suffix) {
        EXPECT_EQ("CGT", Suffix(KMer{"ACGT"}, 1));
        EXPECT_EQ("CGTA", Suffix(KMer{"CGTA"}, 0));
        EXPECT_EQ("", Suffix(KMer{"CC"}, 2));
    }

    TEST(NucleotideToIntTest, NucleotideToInt) {
        struct TestCase {
            char nucleotide;
            int wantResult;
            bool wantSuccess;
        };
        std::vector<TestCase> tests = {
                {'A', 0, true},
                {'C', 1, true},
                {'G', 2, true},
                {'T', 3, true},
                {'B', 0, false},
        };

        for (auto t : tests) {
            int gotResult = 0;
            bool gotSuccess = true;
            try {
                gotResult = NucleotideToInt(t.nucleotide);
            } catch(std::invalid_argument) {
                gotSuccess = false;
            }
            EXPECT_EQ(t.wantSuccess, gotSuccess);
            if (t.wantSuccess) EXPECT_EQ(t.wantResult, gotResult);
        }
    }

    TEST(SuperstringFromPathTest, SuperstringFromPath) {
        struct TestCase {
            std::vector<OverlapEdge> path;
            std::vector<KMer> kMers;
            int k;
            KMerSet wantResult;
        };
        std::vector<TestCase> tests = {
                {
                    std::vector<OverlapEdge>{OverlapEdge{1, 0, 2},OverlapEdge{0, 2, 1}},
                    std::vector<KMer>{KMer{"ACG"}, KMer{"TAC"}, KMer{"GGC"}},
                    3,
                    KMerSet{"TACGGC", std::vector<bool> {1, 1, 0, 1, 0, 0}, 3}
                },
        };

        for (auto t : tests) {
            KMerSet got = SuperstringFromPath(t.path, t.kMers, t.k);
            EXPECT_EQ(t.wantResult.superstring, got.superstring);
            EXPECT_EQ(t.wantResult.k, got.k);
            EXPECT_EQ(t.wantResult.mask.size(), got.mask.size());
            for (int i = 0; i < t.wantResult.mask.size(); ++i) {
                EXPECT_EQ(t.wantResult.mask[i], got.mask[i]);
            }
        }
    }

    TEST(ACAutomatonTest, ConstructTrie) {
        ACAutomaton a;
        std::vector<ACState> wantStates = {
                ACState{std::vector<int>{0, 1}, 1, 0, 0, 0, 0, 0, 0},
                ACState{std::vector<int>{0, 1}, -1, 2, -1, -1, 0, 1, 1},
                ACState{std::vector<int>{0, 1}, -1, -1, 4, 3, 0, 2, 2},
                ACState{std::vector<int>{0}, -1, -1, -1, -1, 0, 3, 3},
                ACState{std::vector<int>{1}, -1, -1, -1, -1, 0, 3, 4},
        };
        std::vector<int> wantEndStateIndices = {3, 4};

        a.ConstructTrie(std::vector<KMer>{KMer{"ACT"}, KMer{"ACG"}});

        ASSERT_EQ(wantStates.size(), a.states.size());
        for (int i = 0; i < wantStates.size(); ++i) {
            EXPECT_EQ(wantStates[i].backwardEdge, a.states[i].backwardEdge);
            for (int j = 0; j < 4; ++j) {
                EXPECT_EQ(wantStates[i].forwardEdges[j], a.states[i].forwardEdges[j]);
            }
            EXPECT_EQ(wantStates[i].depth, a.states[i].depth);
            EXPECT_EQ(wantStates[i].id, a.states[i].id);
            EXPECT_EQ(wantStates[i].supporters, a.states[i].supporters);
        }

        ASSERT_EQ(wantEndStateIndices, a.endStateIndices);
    }

    TEST(ACAutomatonTest, ConstructBackwardEdges) {
        ACAutomaton a;
        // Trie representing ["ACT", "ACA"].
        a.states = {
                ACState{std::vector<int>{0, 1}, 1, 0, 0, 0, 0, 0, 0},
                ACState{std::vector<int>{0, 1}, -1, 2, -1, -1, 0, 1, 1},
                ACState{std::vector<int>{0, 1}, 4, -1, -1, 3, 0, 2, 2},
                ACState{std::vector<int>{0}, -1, -1, -1, -1, 0, 3, 3},
                ACState{std::vector<int>{1}, -1, -1, -1, -1, 0, 3, 4},
        };
        std::vector<ACState> wantStates = {
                ACState{std::vector<int>{0, 1}, 1, 0, 0, 0, 0, 0, 0},
                ACState{std::vector<int>{0, 1}, -1, 2, -1, -1, 0, 1, 1},
                ACState{std::vector<int>{0, 1}, 4, -1, -1, 3, 0, 2, 2},
                ACState{std::vector<int>{0}, -1, -1, -1, -1, 0, 3, 3},
                ACState{std::vector<int>{1}, -1, -1, -1, -1, 1, 3, 4},
        };
        std::vector<int> wantReversedOrdering = {3,4,2,1,0};

        a.ConstructBackwardEdges();

        for (int i = 0; i < wantStates.size(); ++i) {
            EXPECT_EQ(wantStates[i].backwardEdge, a.states[i].backwardEdge);
        }
        EXPECT_EQ(wantReversedOrdering, a.reversedOrdering);
    }

    TEST(GreedyTest, Greedy) {
        std::vector<std::pair<KMerSet, std::vector<KMer>>> tests = {
                {KMerSet{"TACGT",  std::vector<bool> {1, 1, 1, 0, 0}, 3 }, {KMer{"CGT"}, KMer{"TAC"}, KMer{"ACG"}}},
                {KMerSet{"ACGTTT",  std::vector<bool> {1, 1, 0, 1, 0, 0}, 3 }, {KMer{"CGT"}, KMer{"TTT"}, KMer{"ACG"}}},
                {KMerSet{"TACTT",  std::vector<bool> {1, 1, 0, 0, 0}, 4 }, {KMer{"TACT"}, KMer{"ACTT"}}}
        };

        for (auto t : tests) {
            KMerSet got = Greedy(t.second);
            EXPECT_EQ(t.first.superstring, got.superstring);
            EXPECT_EQ(t.first.k, got.k);
            EXPECT_EQ(t.first.mask.size(), got.mask.size());
            for (int i = 0; i < t.first.mask.size(); ++i) {
                EXPECT_EQ(t.first.mask[i], got.mask[i]);
            }
        }
    }
}

int main(int argc, char **argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}