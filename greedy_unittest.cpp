#include "greedy.cpp"

#include "gtest/gtest.h"

namespace {
    TEST(OverlapLengthTest, Trivial) {
        EXPECT_EQ(2, OverlapLength(KMer{"ACGT"}, KMer{"GTAT"}));
        EXPECT_EQ(3, OverlapLength(KMer{"ACCGT"}, KMer{"CGT"}));
        EXPECT_EQ(3, OverlapLength(KMer{"CGT"}, KMer{"CGT"}));
        EXPECT_EQ(3, OverlapLength(KMer{"ACCCTCCCA"}, KMer{"CCATAAAAAACCGGTTAA"}));
    }

    TEST(OverlapLengthTest, Zero) {
        EXPECT_EQ(0, OverlapLength(KMer{"AAA"}, KMer{"TTT"}));
        EXPECT_EQ(0, OverlapLength(KMer{"AAC"}, KMer{"TAA"}));
    }

    TEST(OverlapLengthTest, MultipleOccurrences) {
        EXPECT_EQ(4, OverlapLength(KMer{"TACAC"}, KMer{"ACACGA"}));
        EXPECT_EQ(2, OverlapLength(KMer{"GAA"}, KMer{"AAT"}));
    }

    TEST(Union, Union) {
        EXPECT_EQ("AAAG", Union(KMer{"AAA"}, KMer{"AG"}, 1).value);
        EXPECT_EQ("AAA", Union(KMer{"AAA"}, KMer{"AA"}, 2).value);
        EXPECT_EQ("AAACTG", Union(KMer{"AAA"}, KMer{"CTG"}, 0).value);
    }

    TEST(Greedy, Greedy) {
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