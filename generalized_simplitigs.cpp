#include "models.h"

#include <string>
#include <unordered_set>
#include <vector>
#include <deque>

const char letters[4] {'A', 'C', 'G', 'T'};

/// Convert the encoded KMer representation to string.
std::string NumberToKMer(long long encoded, int length) {
    std::string ret = "";
    for (int i = 0; i < length; ++i) {
        // The last two bits correspond to one nucleotide.
        ret = letters[encoded & 3] + ret;
        // Move to the next letter.
        encoded >>= 2;
    }
    return ret;
}

/// Find the left/right extension to the provided simplitig from the kMers.
/// This extension has k-d overlap with the given simplitig and is left or right based on the fromRight argument.
/// Return the extension - that is the d chars extending the simplitig - and the extending kMer.
std::pair<std::string, std::string> Extension(std::string &simplitig, std::unordered_set<std::string> &kMers, int k, int d, bool fromRight) {
    // Try each of the {A, C, G, T}^d possible extensions of length d.
    for (long long i = 0; i < (1 << (2*d)); ++i) {
        std::string ext = NumberToKMer(i, d);
        std::string next;
        if (fromRight) next = simplitig.substr(simplitig.length() - (k - d),k - d) + ext;
        else next = ext + simplitig.substr(0, k - d);
        if (kMers.count(next) > 0) {
            return {ext, next};
        }
    }
    return {"", ""};
}

/// Find the next generalized simplitig.
/// Update the provided superstring and the mask.
/// Also remove the used k-mers from kMers.
void NextGeneralizedSimplitig(std::unordered_set<std::string> &kMers, std::string &superstring, std::vector<bool> &mask, int k, int d_max) {
    std::string simplitig = *kMers.begin();
    kMers.erase(simplitig);
    std::deque<bool> simplitigMask {1};
    int d_l = 1, d_r = 1;
    while (d_l <= d_max || d_r <= d_max) {
        if (d_r <= d_l) {
            auto extension = Extension(simplitig, kMers, k, d_r, true);
            std::string ext = extension.first;
            if (ext == "") {
                // No right extension found.
                ++d_r;
            } else {
                // Extend the simplitig to the right.
                kMers.erase(extension.second);
                simplitig += ext;
                for (int i = 0; i < d_r - 1; ++i) simplitigMask.push_back(0);
                simplitigMask.push_back(1);
                d_r = 1;
            }
        } else {
            auto extension = Extension(simplitig, kMers, k, d_l, false);
            std::string ext = extension.first;
            if (ext == "") {
                // No left extension found.
                ++d_l;
            } else {
                // Extend the simplitig to the left.
                kMers.erase(extension.second);
                simplitig = ext + simplitig;
                for (int i = 0; i < d_l - 1; ++i) simplitigMask.push_front(0);
                simplitigMask.push_front(1);
                d_l = 1;
            }
        }
    }
    superstring += simplitig;
    for (auto x : simplitigMask) mask.push_back(x);
    for (int i = 0; i < k - 1; ++i) mask.push_back(0);
}

/// Compute the generalized simplitigs greedily.
/// This runs in O(n d_max ^ k), where n is the number of k-mers, but for practical uses is fast for small d_max.
KMerSet GreedyGeneralizedSimplitigs(std::vector<KMer> kMers, int k, int d_max) {
    std::string superstring = "";
    std::vector<bool> mask;

    std::unordered_set<std::string> remainingKMers;
    for (auto &&kMer : kMers) remainingKMers.insert(kMer.value);
    while(!remainingKMers.empty()) NextGeneralizedSimplitig(remainingKMers, superstring, mask, k, d_max);
    return KMerSet{superstring, mask, k};
}
