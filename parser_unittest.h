#pragma once
#include "parser.h"

#include <algorithm>
#include <vector>
#include <string>

#include "gtest/gtest.h"

namespace {
    TEST(ReadFastaTest, ReadFasta) {
        std::string path = (std::string) get_current_dir_name();
        path += "/tests/test.fa";
        std::vector<FastaRecord> wantResult = {
                FastaRecord{">1", "ACCCGAAC"},
                FastaRecord{">2", "CGTANATGC"},
                FastaRecord{">3 with description", "ACCCGTTTAACG"},
                FastaRecord{">4", "A"},
        };

        auto gotResult = ReadFasta(path);

        EXPECT_EQ(wantResult.size(), gotResult.size());
        for (size_t i = 0; i < wantResult.size(); ++i) {
            EXPECT_EQ(wantResult[i].name, gotResult[i].name);
            EXPECT_EQ(wantResult[i].sequence, gotResult[i].sequence);
        }
    }

    TEST(ReadKMersTest, ReadKMers) {
        std::string path = (std::string) get_current_dir_name();
        path += "/tests/test.fa";
        std::vector<FastaRecord> wantResult = {
                FastaRecord{">1", "ACCCGAAC"},
                FastaRecord{">2", "CGTANATGC"},
                FastaRecord{">3 with description", "ACCCGTTTAACG"},
                FastaRecord{">4", "A"},
        };
        struct TestCase {
            int k;
            bool complements;
            size_t wantResultSize;
        };
        std::vector<TestCase> tests = {
                {10, true, 3},
                {5, true, 11},
                {5, false, 11},
                {2, true, 9},
                {2, false, 11},
        };

        for (auto &t: tests) {
            auto gotResult = ReadKMers(path, t.k, t.complements);
            EXPECT_EQ(t.wantResultSize, gotResult.size());
        }

    }

    TEST(AddKMersFromSequenceTest, String) {
        struct TestCase {
            std::string data;
            std::unordered_set<std::string> initialKMers;
            int k;
            std::unordered_set<std::string> wantResult;
        };
        std::vector<TestCase> tests = {
                {"AAAA", {}, 2, {{"AA"}}},
                {"ACGTA", {}, 3, {{"ACG"}, {"CGT"}, {"GTA"}}},
                {"GAAAAGTTTAAAAAGAC", {},  4,
                 {{"AAAA"}, {"AAAG"}, {"AAGA"}, {"AAGT"},
                                   {"AGAC"}, {"AGTT"}, {"GAAA"}, {"GTTT"}, {"TAAA"},
                                   {"TTAA"}, {"TTTA"}}},
                {"ACGTA", {"ACG", "CCC"}, 3, {{"CCC"}, {"ACG"}, {"CGT"}, {"GTA"}}},
                {"acgTa", {}, 3, {{"ACG"}, {"CGT"}, {"GTA"}}},
                {"ACGNTTA", {}, 3, {{"ACG"}, {"TTA"}}},
                {"ACGNMTTA", {}, 3, {{"ACG"}, {"TTA"}}},
        };

        for (auto t: tests) {
            AddKMersFromSequence(t.initialKMers, t.data, t.k);
            for (auto &&kMer : t.initialKMers) EXPECT_EQ(1, t.wantResult.count(kMer));
            for (auto &&kMer : t.wantResult) EXPECT_EQ(1, t.initialKMers.count(kMer));
        }
    }

    TEST(AddKMersFromSequenceTest, Number) {
        struct TestCase {
            std::string data;
            std::unordered_set<int64_t> initialKMers;
            int k;
            bool complements;
            size_t wantResultSize;
        };
        std::vector<TestCase> tests = {
                {"AAAA", {}, 2, true, 1},
                {"ACGTA", {}, 3, true, 2},
                {"ACGTA", {}, 3, false, 3},
                {"acgTa", {}, 3, false, 3},
                {"ACGNMTTA", {}, 3, true, 2},
                {"ACCCGTTTAACG", {}, 10, true, 3},
                {"ACC", {}, 10, true, 0},

        };
        for (auto t: tests) {
            AddKMersFromSequence(t.initialKMers, t.data, t.k, t.complements);
            EXPECT_EQ(t.wantResultSize, t.initialKMers.size());
        }
    }

    TEST(FilterKMersWithComplementTest, FilterKMersWithComplement) {
        struct TestCase {
            std::unordered_set<std::string> kMers;
            size_t wantResultSize;
        };
        std::vector<TestCase> tests = {
                {{"AA", "TT", "AC", "GT"}, 2},
                {{"AT", }, 1},
                {{"GT", "CG"}, 2},
        };

        for (auto t: tests) {
            auto gotResult = FilterKMersWithComplement(t.kMers);
            EXPECT_EQ(t.wantResultSize, gotResult.size());
        }
    }

    TEST(ConstructKMersTest, ConstructKMers) {
        struct TestCase {
            std::vector<FastaRecord> data;
            int k;
            std::vector<KMer> wantResult;
            bool complements;
        };
        std::vector<TestCase> tests = {
                {{FastaRecord{"", "AAAA"}}, 2, std::vector<KMer>{KMer{"AA"}}, false},
                {{FastaRecord{"", "GAAAAGTTTAAAAAGAC"}}, 4,
                        std::vector<KMer>{KMer{"AAAA"}, KMer{"AAAG"}, KMer{"AAGA"}, KMer{"AAGT"},
                                          KMer{"AGAC"}, KMer{"AGTT"}, KMer{"GAAA"}, KMer{"GTTT"}, KMer{"TAAA"},
                                          KMer{"TTAA"}, KMer{"TTTA"}}, false},
                {{FastaRecord{"", "AAAA"}, FastaRecord{"", "CCC"}}, 2, std::vector<KMer>{KMer{"AA"}, KMer{"CC"}}, false},
                {{FastaRecord{"", "AATT"}}, 2, std::vector<KMer>{KMer{"AT"}, KMer{"TT"}}, true},
        };

        for (auto t: tests) {
            auto gotResult = ConstructKMers(t.data, t.k, t.complements);
            // ConstructKMers does not return the kMers in any particular order.
            std::sort(gotResult.begin(), gotResult.end(), [](const KMer &a, const KMer &b) {return a.value < b.value;});
            EXPECT_EQ(t.wantResult.size(), gotResult.size());
            for (size_t i = 0; i < t.wantResult.size(); ++i) {
                EXPECT_EQ(t.wantResult[i].value, gotResult[i].value);
            }
        }
    }
}
