#pragma once
#include "../src/masks.h"

#include "kmer_types.h"

#include <algorithm>
#include <vector>
#include <string>
#include <filesystem>

#include "gtest/gtest.h"

namespace {

// Retrieving current path on Windows does not work as on linux
// therefore the following unittest is linux-specific.
#ifdef __unix__

    TEST(Masks, OptimizeOnes) {
        std::string path = std::filesystem::current_path();
        path += "/tests/testdata/masktest.fa";

        struct TestCase {
            std::vector<kmer_t> kMers;
            int k;
            bool complements;
            bool minimize;
            std::string wantResult;
        };
        std::vector<TestCase> tests = {
                {
                    {KMerToNumber({"ACG"}), KMerToNumber({"CGT"}), KMerToNumber({"TAA"})},
                    3, false, true,
                    ">superstring reoptimized=minone\nACgTaacgt\n"
                },
                {
                    {KMerToNumber({"ACG"}), KMerToNumber({"CGT"}), KMerToNumber({"TAA"})},
                    3, false, false,
                    ">superstring reoptimized=maxone\nACgTaACgt\n"
                },
                {
                    {KMerToNumber({"ACG"}), KMerToNumber({"TAA"})},
                    3, true, true,
                    ">superstring reoptimized=minone\nAcgTaacgt\n"
                },
                {
                    {KMerToNumber({"ACG"}), KMerToNumber({"TAA"})},
                    3, true, false,
                    ">superstring reoptimized=maxone\nACgTaACgt\n"
                },
        };

        for (auto &t : tests) {
            std::stringstream of;
            auto masked_superstring = ReadMaskedSuperstring(path);
            auto kMersDict = wrapper.kh_init_set();
            int ret;
            for (auto &kMer : t.kMers) wrapper.kh_put_to_set(kMersDict, kMer, &ret);

            OptimizeOnes(masked_superstring, of, kMersDict, wrapper, kmer_t (0), t.k, t.complements, t.minimize);
            kseq_destroy(masked_superstring);

            EXPECT_EQ(t.wantResult, of.str());
        }
    }

    TEST(Mask, OptimizeRuns) {
        std::string path = std::filesystem::current_path();

        struct TestCase {
            std::vector<kmer_t> kMers;
            int k;
            bool complements;
            bool approximate;
            std::string relativePath;
            std::string wantResult;
        };
        std::vector<TestCase> tests = {
                {
                        {
                            KMerToNumber({"AC"}),
                            KMerToNumber({"CG"}),
                            KMerToNumber({"GT"}),
                            KMerToNumber({"TT"}),
                            KMerToNumber({"AT"})
                        },
                        2, false, false,
                        "/tests/testdata/runstest.fa",
                        ">superstring reoptimized=minrun\nacgcgttACGtATt\n"
                },
                {
                        {
                                KMerToNumber({"AC"}),
                                KMerToNumber({"CG"}),
                                KMerToNumber({"GT"}),
                                KMerToNumber({"TT"}),
                                KMerToNumber({"AT"})
                        },
                        2, false, true,
                        "/tests/testdata/runstest.fa",
                        ">superstring reoptimized=approxminrun\nACgCGTtACGtATt\n"
                },
                {
                        {
                                KMerToNumber({"ACT"}),
                                KMerToNumber({"CTA"}),
                                KMerToNumber({"GTA"})
                        },
                        3, true, false,
                        "/tests/testdata/runstest2.fa",
                        ">superstring reoptimized=minrun\nACTAGta\n"
                },
        };

        for (auto &t : tests) {
            std::stringstream of;
            auto totalPath = path + t.relativePath;
            auto masked_superstring = ReadMaskedSuperstring(totalPath);
            auto kMersDict = wrapper.kh_init_set();
            int ret;
            for (auto &kMer : t.kMers) wrapper.kh_put_to_set(kMersDict, kMer, &ret);

            OptimizeRuns(wrapper, kmer_t (0), masked_superstring, kMersDict, of, t.k, t.complements, t.approximate);

            EXPECT_EQ(t.wantResult, of.str());
        }
    }
#endif
}
