#include "models.h"
#include <string>
#include <unordered_set>
#include <vector>
#include <deque>

const char letters[4] {'A', 'C', 'G', 'T'};

std::string NumberToKMer(long long encoded, int length) {
    std::string ret = "";
    for (int i = 0; i < length; ++i) {
        ret = letters[encoded & 3] + ret;
        encoded >>= 2;
    }
    return ret;
}

std::pair<std::string, std::string> Extension(std::string &simplitig, std::unordered_set<std::string> &kMers, int k, int d, bool fromRight) {
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
                ++d_r;
            } else {
                kMers.erase(extension.second);
                simplitig += ext;
                for (int i = 0; i < d_r - 1; ++i) simplitigMask.push_back(0);
                simplitigMask.push_back(1);
            }
        } else {
            auto extension = Extension(simplitig, kMers, k, d_l, false);
            std::string ext = extension.first;
            if (ext == "") {
                ++d_l;
            } else {
                kMers.erase(extension.second);
                simplitig = ext + simplitig;
                for (int i = 0; i < d_r - 1; ++i) simplitigMask.push_front(0);
                simplitigMask.push_front(1);
            }
        }
    }
    superstring += simplitig;
    for (auto x : simplitigMask) mask.push_back(x);
    for (int i = 0; i < k - 1; ++i) mask.push_back(0);
}

KMerSet GreedyGeneralizedSimplitigs(std::vector<KMer> kMers, int k, int d_max) {
    std::string superstring = "";
    std::vector<bool> mask;

    std::unordered_set<std::string> remainingKMers;
    for (auto &&kMer : kMers) remainingKMers.insert(kMer.value);
    while(!remainingKMers.empty()) NextGeneralizedSimplitig(remainingKMers, superstring, mask, k, d_max);
    return KMerSet{superstring, mask, k};
}

