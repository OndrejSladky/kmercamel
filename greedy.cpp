#include "models.h"

#include <vector>
#include <iostream>
#include <unordered_map>
#include <list>

#include "greedy_ac.cpp"

long long KMerToNumber(KMer &kMer) {
    long long ret = 0;
    for (int i = kMer.length() - 1; i >= 0; --i) {
        ret <<= 2;
        ret |= NucleotideToInt(kMer.value[i]);
    }
    return ret;
}

long long BitPrefix(long long kMer, int k, int d) {
    long long mask = -1LL ^ ((1 << ((k - d) << 1)) - 1);
    return (kMer & mask) >> ((k - d) << 1);
}

long long BitSuffix(long long kMer, int d) {
    return kMer & ((1 << (d << 1)) - 1);
}

/// Greedily find the approximate hamiltonian path with longest overlaps.
std::vector<OverlapEdge> OverlapHamiltonianPath (std::vector<KMer> &input, int k) {
    std::vector<long long> kMers(input.size());
    std::vector<OverlapEdge> hamiltonianPath;
    std::vector<bool> suffixForbidden(kMers.size(), false);
    std::vector<bool> prefixForbidden(kMers.size(), false);
    std::vector<int> first(kMers.size());
    std::vector<int> last(kMers.size());
    for (int i = 0; i < kMers.size(); ++i) {
        first[i] = last[i] = i;
        kMers[i] = KMerToNumber(input[i]);
    }
    for (int d = k - 1; d >= 0; --d) {
        std::unordered_map<long long, std::list<int>> prefixes;
        for (int i = 0 ; i < kMers.size(); ++i) if(!prefixForbidden[i]) {
            long long prefix = BitPrefix(kMers[i], k, d);
            if (prefixes.count(prefix) == 0) prefixes[prefix] = std::list<int>();
            prefixes[prefix].push_back(i);
        }
        for (int i = 0 ; i < kMers.size(); ++i) if(!suffixForbidden[i]) {
            long long suffix = BitSuffix(kMers[i], d);
            if (prefixes.count(suffix) == 0 || prefixes[suffix].size() == 0) continue;
            auto j = prefixes[suffix].begin();
            if (first[*j] == i) {
                if (prefixes[suffix].size() == 1) continue;
                j++;
            }
            hamiltonianPath.push_back(OverlapEdge{*j, i, d});
            suffixForbidden[i] = true;
            prefixForbidden[*j] = true;
            first[last[i]] = first[*j];
            last[first[*j]] = last[i];
            prefixes[suffix].erase(j);
        }
    }
    return hamiltonianPath;
}

/// Get the approximated shortest superstring of the given k-mers using the GREEDY algorithm.
/// This runs in O(n k), where n is the number of k-mers.
KMerSet Greedy(std::vector<KMer> &kMers) {
    if (kMers.size() == 0) {
        throw new std::invalid_argument("input cannot be empty");
    }
    int k = kMers[0].length();
    auto hamiltonianPath = OverlapHamiltonianPath(kMers, k);
    return SuperstringFromPath(hamiltonianPath, kMers, k);
}
