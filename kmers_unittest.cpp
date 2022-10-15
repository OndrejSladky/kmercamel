#include "kmers.cpp"

#include <algorithm>
#include <vector>
#include <string>

#include "gtest/gtest.h"

namespace {
    TEST(ConstructKMersTest, ConstructKMers) {
        struct TestCase {
            std::string data;
            int k;
            std::vector<KMer> wantResult;
        };
        std::vector<TestCase> tests = {
                {"AAAA", 2, std::vector<KMer>{KMer{"AA"}}},
                {"ACGTA", 3, std::vector<KMer>{KMer{"ACG"}, KMer{"CGT"}, KMer{"GTA"}}},
                {"GAAAAGTTTAAAAAGAC", 4,
                 std::vector<KMer>{KMer{"AAAA"}, KMer{"AAAG"}, KMer{"AAGA"}, KMer{"AAGT"},
                                   KMer{"AGAC"}, KMer{"AGTT"}, KMer{"GAAA"}, KMer{"GTTT"}, KMer{"TAAA"},
                                   KMer{"TTAA"}, KMer{"TTTA"}}},
        };

        for (auto t: tests) {
            auto gotResult = ConstructKMers(t.data, t.k);
            // ConstructKMers does not return the kMers in any particular order.
            std::sort(gotResult.begin(), gotResult.end(), [](const KMer &a, const KMer &b) {return a.value < b.value;});
            EXPECT_EQ(t.wantResult.size(), gotResult.size());
            for (int i = 0; i < t.wantResult.size(); ++i) {
                EXPECT_EQ(t.wantResult[i].value, gotResult[i].value);
            }
        }
    }
}
