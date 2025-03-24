#pragma once
#include "../src/simplitigs.h"


#include "kmer_types.h"

#include "gtest/gtest.h"
typedef unsigned char byte;
namespace {
    TEST(Simplitigs, SimplitigReverseComplement) {
        struct TestCase {
            simplitig_t input;
            simplitig_t want_result;
        };
        std::vector<TestCase> tests = {
                {
                        {0, 0, 1, 0, 1, 0, 0, 1}, // AGGC
                        {1, 0, 0, 1, 0, 1, 1, 1}, // GCCT
                },
                {
                        {0, 0, 1, 1, 1, 0}, // ATG
                        {0, 1, 0, 0, 1, 1}, // CAT
                },
        };

        for (auto t : tests) {
            simplitig_t got_result = simplitig_reverse_complement(t.input);

            EXPECT_EQ(t.want_result, got_result);
        }
    }
    
    TEST(Simplitigs, SimplitigAtIndex) {
        struct TestCase {
            simplitig_t simplitig;
            size_t index;
            char want_result;
        };
        std::vector<TestCase> tests = {
                {
                        {0, 0, 1, 0, 1, 0, 0, 1}, // AGGC
                        0, 'A'
                },
                {
                        {0, 0, 1, 0, 1, 0, 0, 1}, // AGGC
                        1, 'G'
                },
                {
                        {0, 0, 1, 0, 1, 0, 0, 1}, // AGGC
                        2, 'G'
                },
                {
                        {0, 0, 1, 0, 1, 0, 0, 1}, // AGGC
                        3, 'C'
                },
                {
                        {0, 0, 1, 1, 1, 0}, // ATG
                        1, 'T'
                },
        };

        for (auto t : tests) {
            char got_result = simplitig_at_index(t.simplitig, t.index);

            EXPECT_EQ(t.want_result, got_result);
        }
    }

    TEST(Simplitigs, KmersInSimplitig) {
        struct TestCase {
            simplitig_t simplitig;
            int k;
            size_t want_result;
        };
        std::vector<TestCase> tests = {
                {
                        {0, 0, 1, 0, 1, 0, 0, 1}, // AGGC
                        1, 4
                },
                {
                        {0, 0, 1, 0, 1, 0, 0, 1}, // AGGC
                        2, 3
                },
                {
                        {0, 0, 1, 0, 1, 0, 0, 1}, // AGGC
                        3, 2
                },
                {
                        {0, 0, 1, 0, 1, 0, 0, 1}, // AGGC
                        4, 1
                },
                {
                        {0, 0, 1, 1, 1, 0}, // ATG
                        2, 2
                },
        };

        for (auto t : tests) {
            size_t got_result = kmers_in_simplitig(t.simplitig, t.k);

            EXPECT_EQ(t.want_result, got_result);
        }
    }

    TEST(Simplitigs, KmerAtSimplitigIndex) {
        struct TestCase {
            simplitig_t simplitig;
            int k;
            size_t index;
            kmer_t want_result;
        };
        std::vector<TestCase> tests = {
                {
                        {0, 0, 1, 0, 1, 0, 0, 1}, // AGGC
                        2, 0, 0b0010,
                },
                {
                        {0, 0, 1, 0, 1, 0, 0, 1}, // AGGC
                        2, 1, 0b1010
                },
                {
                        {0, 0, 1, 0, 1, 0, 0, 1}, // AGGC
                        2, 2, 0b1001
                },
                {
                        {0, 0, 1, 0, 1, 0, 0, 1}, // AGGC
                        3, 1, 0b101001
                },
                {
                        {0, 0, 1, 1, 1, 0}, // ATG
                        3, 0, 0b001110
                },
        };

        for (auto t : tests) {
            kmer_t got_result = kmer_at_simplitig_index(kmer_t(0), t.simplitig, t.k, t.index);

            EXPECT_EQ(t.want_result, got_result);
        }
    }

    TEST(Local, NextSimplitig) {
        struct TestCase {
            std::vector<kmer_t> kmers;
            int k;
            bool complements;
            simplitig_t want_simplitig;
            // Sorted k-mers.
            std::vector<kmer_t> want_kmers;
        };
        std::vector<TestCase> tests = {
                // {ACAA, AACA}
                {{0b00010000,  0b00000100}, 4, false,
                 simplitig_from_string("AACAA"), {},},
                // {ACAA, ATTT, TGTT, AAAT, TTGT, AACA}
                {{0b00010000, 0b00111111, 0b11101111, 0b00000011, 0b11111011, 0b00000100}, 4, true,
                 simplitig_from_string("AACAA"), {0b00000011, 0b00111111},},
        };

        for (auto &&t: tests) {
            std::stringstream of;
            auto kmers = wrapper.kh_init_set();
            int ret;
            for (auto &&kmer : t.kmers) wrapper.kh_put_to_set(kmers, kmer, &ret);

            simplitig_t got_simplitig;
            if (t.complements) got_simplitig = next_simplitig<true>(kmers, wrapper, t.kmers.front(), t.k);
            else got_simplitig = next_simplitig<false>(kmers, wrapper, t.kmers.front(), t.k);
            
            auto remaining_kmers = kMersToVec(kmers, kmer_t(0));
            std::sort(remaining_kmers.begin(), remaining_kmers.end());

            EXPECT_EQ(t.want_kmers, remaining_kmers);
            // Check that a valid superstring was returned.
            bool valid = t.want_simplitig == got_simplitig;
            if (t.complements) valid |= simplitig_reverse_complement(t.want_simplitig) == got_simplitig;
            if (!valid) {
                EXPECT_EQ(t.want_simplitig, got_simplitig);
            }
        }
    }

    TEST(Simplitigs, GetSimplitigs) {
        struct TestCase {
            std::vector<kmer_t> kmers;
            int k;
            bool complements;
            std::vector<simplitig_t> want_simplitigs;
        };
        std::vector<TestCase> tests = {
                {{KMerToNumber(KMer{"GCT"}), KMerToNumber(KMer{"TAA"}), KMerToNumber(KMer{"AAA"})}, 3, false, {simplitig_from_string("GCT"), simplitig_from_string("TAAA")}},
                {{KMerToNumber(KMer{"TAA"}), KMerToNumber(KMer{"AAA"}), KMerToNumber(KMer{"GCT"})}, 3, false, {simplitig_from_string("GCT"), simplitig_from_string("TAAA")}},
                {{KMerToNumber(KMer{"TTTCTTTTTTTTTTTTTTTTTTTTTTTTTTG"}), KMerToNumber(KMer{"TTCTTTTTTTTTTTTTTTTTTTTTTTTTTGA"})}, 31, false,
                 {simplitig_from_string("TTtcttttttttttttttttttttttttttga")}},
        };

        for (auto t: tests) {
            std::stringstream of;
            auto kmers = wrapper.kh_init_set();
            int ret;
            for (auto &&kmer : t.kmers) wrapper.kh_put_to_set(kmers, kmer, &ret);

            auto got_result = get_simplitigs(kmers, wrapper, kmer_t(0), t.k, t.complements);
            sort(got_result.begin(), got_result.end(), [](const simplitig_t &a, const simplitig_t &b) {return a.size() < b.size();});

            EXPECT_EQ(t.want_simplitigs, got_result);
        }
    }
}
