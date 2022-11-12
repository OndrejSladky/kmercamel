#include "models.h"

#include <vector>
#include <iostream>
#include <unordered_map>
#include <list>
#include <cstdint>

#include "greedy_ac.cpp"
#include "kmers.cpp"

/// Greedily find the approximate hamiltonian path with longest overlaps.
/// k is the size of one k-mer and n is the number of distinct k-mers.
/// If complements are provided, treat k-mer and its complement as identical.
/// If this is the case, k-mers are expected not to contain both k-mer and its complement.
std::vector<OverlapEdge> OverlapHamiltonianPath (std::vector<KMer> &input, int k, bool complements) {
    // If complements, kMers should be twice the size.
    std::vector<int64_t> kMers((1 + complements) * input.size());
    std::vector<OverlapEdge> hamiltonianPath;
    std::vector<bool> suffixForbidden(kMers.size(), false);
    std::vector<bool> prefixForbidden(kMers.size(), false);
    std::vector<size_t> first(kMers.size());
    std::vector<size_t> last(kMers.size());
    for (size_t i = 0; i < input.size(); ++i) {
        first[i] = last[i] = i;
        kMers[i] = KMerToNumber(input[i]);
    }
    // Fill in the reverse complements. K-mer and its reverse complement are at indices i and i + input.size.
    for (size_t i = input.size(); i < kMers.size(); ++i) {
        first[i] = last[i] = i;
        kMers[i] = ReverseComplement(KMerToNumber(input[i - input.size()]), k);
    }
    for (int d = k - 1; d >= 0; --d) {
        std::unordered_map<int64_t, std::list<size_t>> prefixes;
        for (size_t i = 0 ; i < kMers.size(); ++i) if(!prefixForbidden[i]) {
            int64_t prefix = BitPrefix(kMers[i], k, d);
            if (prefixes.count(prefix) == 0) prefixes[prefix] = std::list<size_t>();
            prefixes[prefix].push_back(i);
        }
        for (size_t i = 0 ; i < kMers.size(); ++i) if(!suffixForbidden[i]) {
            int64_t suffix = BitSuffix(kMers[i], d);
            if (prefixes.count(suffix) == 0 || prefixes[suffix].empty()) continue;
            auto j = prefixes[suffix].begin();
            // If the path forms a cycle, or is between k-mer and its reverse complement, skip this path.
            while (j != prefixes[suffix].end() && (first[i] == *j || (i % input.size()) == (*j % input.size()) || prefixForbidden[*j])) {
                auto new_j = j;
                new_j++;
                if (prefixForbidden[*j]) prefixes[suffix].erase(j);
                j = new_j;
            }
            if (j == prefixes[suffix].end()) {
                continue;
            }
            hamiltonianPath.push_back(OverlapEdge{i, *j, d});
            suffixForbidden[i] = true;
            prefixForbidden[*j] = true;
            // Forbid the reverse complements from forming any paths.
            if (complements) {
                size_t i_comp = (i + input.size()) % kMers.size();
                size_t j_comp = (*j + input.size()) % kMers.size();
                suffixForbidden[i_comp] = true;
                suffixForbidden[j_comp] = true;
                prefixForbidden[i_comp] = true;
                prefixForbidden[j_comp] = true;
            }
            first[last[*j]] = first[i];
            last[first[i]] = last[*j];
            prefixes[suffix].erase(j);
        }
    }
    return hamiltonianPath;
}

/// Get the approximated shortest superstring of the given k-mers using the GREEDY algorithm.
/// This runs in O(n k), where n is the number of k-mers.
/// If complements are provided, treat k-mer and its complement as identical.
/// If this is the case, k-mers are expected not to contain both k-mer and its complement.
KMerSet Greedy(std::vector<KMer> kMers, bool complements) {
    if (kMers.empty()) {
        throw std::invalid_argument("input cannot be empty");
    }
    int k = kMers[0].length();
    auto hamiltonianPath = OverlapHamiltonianPath(kMers, k, complements);
    if (complements) {
      size_t n = kMers.size();
      kMers.resize(2 * n);
      for (size_t i = 0; i < n; ++i) {
          kMers[i + n] = ReverseComplement(kMers[i]);
      }
    }
    return SuperstringFromPath(hamiltonianPath, kMers, k);
}
