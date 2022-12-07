#pragma once
#include "models.h"
#include "kmers.h"

#include <string>
#include <unordered_set>
#include <vector>
#include <deque>
#include <cstdint>


/// Find the right extension to the provided last k-mer from the kMers.
/// This extension has k-d overlap with the given simplitig.
/// Return the extension - that is the d chars extending the simplitig - and the extending kMer.
std::pair<int64_t, int64_t> RightExtension(int64_t last, std::unordered_set<int64_t> &kMers, int k, int d) {
    // Try each of the {A, C, G, T}^d possible extensions of length d.
    for (int64_t ext = 0; ext < (1 << (d << 1)); ++ext) {
        int64_t next = BitSuffix(last, k - d) << (d << 1) | ext;
        if (kMers.count(next) > 0) {
            return {ext, next};
        }
    }
    return {-1, -1};
}

/// Find the left extension to the provided first k-mer from the kMers.
/// This extension has k-d overlap with the given simplitig.
/// Return the extension - that is the d chars extending the simplitig - and the extending kMer.
std::pair<int64_t, int64_t> LeftExtension(int64_t first, std::unordered_set<int64_t> &kMers, int k, int d) {
    // Try each of the {A, C, G, T}^d possible extensions of length d.
    for (int64_t ext = 0; ext < (1 << (d << 1)); ++ext) {
        int64_t next = ext << ((k - d) << 1) | BitPrefix(first, k, k - d);
        if (kMers.count(next) > 0) {
            return {ext, next};
        }
    }
    return {-1, -1};
}

/// Find the next generalized simplitig.
/// Update the provided superstring and the mask.
/// Also remove the used k-mers from kMers.
/// If complements are true, it is expected that kMers contain both k-mer and its reverse complement.
void NextGeneralizedSimplitig(std::unordered_set<int64_t> &kMers, std::string &superstring, std::vector<bool> &mask, int k, int d_max, bool complements) {
     // Maintain the first and last k-mer in the simplitig.
    int64_t last = *kMers.begin(), first = last;
    std::string simplitig = NumberToKMer(last, k);
    kMers.erase(last);
    if (complements) kMers.erase(ReverseComplement(last, k));
    std::deque<bool> simplitigMask {1};
    int d_l = 1, d_r = 1;
    while (d_l <= d_max || d_r <= d_max) {
        if (d_r <= d_l) {
            auto extension = RightExtension(last, kMers, k, d_r);
            int64_t ext = extension.first;
            if (ext == -1) {
                // No right extension found.
                ++d_r;
            } else {
                // Extend the simplitig to the right.
                kMers.erase(extension.second);
                if (complements) kMers.erase(ReverseComplement(extension.second, k));
                simplitig += NumberToKMer(ext, d_r);
                for (int i = 0; i < d_r - 1; ++i) simplitigMask.push_back(0);
                simplitigMask.push_back(1);
                d_r = 1;
                last = extension.second;
            }
        } else {
            auto extension = LeftExtension(first, kMers, k, d_l);
            int64_t ext = extension.first;
            if (ext == -1) {
                // No left extension found.
                ++d_l;
            } else {
                // Extend the simplitig to the left.
                kMers.erase(extension.second);
                if (complements) kMers.erase(ReverseComplement(extension.second, k));
                simplitig = NumberToKMer(ext, d_l) + simplitig;
                for (int i = 0; i < d_l - 1; ++i) simplitigMask.push_front(0);
                simplitigMask.push_front(1);
                d_l = 1;
                first = extension.second;
            }
        }
    }
    superstring += simplitig;
    for (auto x : simplitigMask) mask.push_back(x);
    // Fill the remaining zeros of the last k-mer in the simplitig.
    for (int i = 0; i < k - 1; ++i) mask.push_back(0);
}

/// Compute the generalized simplitigs greedily.
/// This runs in O(n d_max ^ k), where n is the number of k-mers, but for practical uses it is fast.
KMerSet GreedyGeneralizedSimplitigs(std::vector<KMer> kMers, int k, int d_max, bool complements) {
    std::string superstring;
    std::vector<bool> mask;

    std::unordered_set<int64_t> remainingKMers;
    for (auto &&kMer : kMers) remainingKMers.insert(KMerToNumber(kMer));
    if (complements) for (auto &&kMer : kMers) remainingKMers.insert(ReverseComplement(KMerToNumber(kMer), k));
    while(!remainingKMers.empty()) NextGeneralizedSimplitig(remainingKMers, superstring, mask, k, d_max, complements);
    return KMerSet{superstring, mask, k};
}
