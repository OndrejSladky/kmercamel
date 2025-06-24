#pragma once

#include "../src/streaming.h"
#include "gtest/gtest.h"
namespace {

// Retrieving current path on Windows does not work as on linux
// therefore the following unittest is linux-specific.
#ifdef __unix__
    TEST(Streaming, Streaming) {
        std::string path = std::filesystem::current_path();
        path += "/tests/testdata/test.fa";
        struct TestCase {
            int k;
            bool complements;
            std::string wantResult;
        };
        std::vector<TestCase> tests = {
                {3, true, "ACCCGAacCGtaATgcTTta"},
                {2, true, "ACcCGAaTaATGc"},
                {1, true, "AC"},
                {1, false, "ACGT"},
        };


        for (auto t : tests) {
            std::stringstream of;

            Streaming(kmer_dict64_t(), kmer64_t(0), path, of, t.k, t.complements);

            EXPECT_EQ(t.wantResult, of.str());
        }
    }
    TEST(Streaming, StreamingFiltered) {
        std::string path = std::filesystem::current_path();
        path += "/tests/testdata/test.fa";
        struct TestCase {
            int k;
            bool complements;
            uint16_t min_frequency;
            std::string wantResult;
        };
        std::vector<TestCase> tests = {
                {3, true, 2, "ACCCGttTaa"},
                {3, false, 2, "ACCCgtAac"},
                {3, true, 3, "AAcg"},
                {4, true, 2, "ACccgAacg"},
                {1, true, 2, "CA"},
                {1, false, 2, "CAGT"},
                {1, false, 3, "CAGT"},
                {1, false, 4, "CAGT"},
                {1, false, 5, "CATG"},
                {1, false, 6, "CA"},
                {1, false, 10, "C"},
        };

        for (auto t : tests) {
            std::stringstream of;

            StreamingFiltered(kmer_dict64_t(), kmer64_t(0), path, of, t.k, t.complements, t.min_frequency);

            EXPECT_EQ(t.wantResult, of.str());
        }
    }
#endif
}
