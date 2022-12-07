#pragma once
#include "models.h"

#include <vector>
#include <iostream>
#include <unordered_map>
#include <list>
#include <cstdint>

#include "kmers.h"

/// Represents oriented edge in the overlap graph.
struct OverlapEdge {
    // Index of the first k-mer.
    size_t firstIndex;
    // Index of the second k-mer.
    size_t secondIndex;
    // Length of the overlap of the two k-mers.
    int overlapLength;
};

/// Greedily find the approximate Hamiltonian path with longest overlaps.
/// k is the size of one k-mer and n is the number of distinct k-mers.
/// If complements are provided, treat k-mer and its complement as identical.
/// If this is the case, k-mers are expected to contain both the k-mer and its reverse complement at positions i and i+size/2.
/// Moreover, if so, the resulting Hamiltonian path contains two superstrings which are reverse complements of one another.
std::vector<OverlapEdge> OverlapHamiltonianPath (std::vector<int64_t> &kMers, int k, bool complements) {
    size_t n = kMers.size() / (1 + complements);
    std::vector<OverlapEdge> hamiltonianPath;
    std::vector<bool> suffixForbidden(kMers.size(), false);
    std::vector<bool> prefixForbidden(kMers.size(), false);
    std::vector<size_t> first(kMers.size());
    std::vector<size_t> last(kMers.size());
    for (size_t i = 0; i < kMers.size(); ++i) {
        first[i] = last[i] = i;
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
            // If the path forms a cycle, or is between k-mer and its reverse complement, or the k-mers complement was already selected skip this path.
            while (j != prefixes[suffix].end() && (first[i]%n == *j%n || first[i]%n == last[*j]%n || prefixForbidden[*j])) {
                auto new_j = j;
                new_j++;
                // If the k-mer is forbidden, remove it to keep the complexity linear.
                if (prefixForbidden[*j]) prefixes[suffix].erase(j);
                j = new_j;
            }
            if (j == prefixes[suffix].end()) {
                continue;
            }
            std::vector<std::pair<size_t,size_t>> new_edges ({{i, *j}});
            // Add also the edge between complementary k-mers in the opposite direction.
            if (complements) new_edges.emplace_back((*j + n) % kMers.size(), (i + n) % kMers.size());
            for (auto [x, y] : new_edges) {
                hamiltonianPath.push_back(OverlapEdge{x, y, d});
                prefixForbidden[y] = true;
                first[last[y]] = first[x];
                last[first[x]] = last[y];
                suffixForbidden[x] = true;
            }
            prefixes[suffix].erase(j);
        }
    }
    return hamiltonianPath;
}

/// Construct the superstring and its mask from the given hamiltonian path in the overlap graph.
/// If reverse complements are considered and the hamiltonian path contains two paths which are reverse complements of one another,
/// return only one of them.
KMerSet SuperstringFromPath(const std::vector<OverlapEdge> &hamiltonianPath, const std::vector<int64_t> &kMers, const int k) {
    std::vector<OverlapEdge> edgeFrom (kMers.size(), OverlapEdge{size_t(-1),size_t(-1), -1});
    std::vector<bool> isStart(kMers.size(), true);
    for (auto edge : hamiltonianPath) {
        edgeFrom[edge.firstIndex] = edge;
    }
    for (auto edge : hamiltonianPath) {
        isStart[edge.secondIndex] = false;
    }

    // Find the vertex in the overlap graph with in-degree 0.
    size_t start = 0;
    for (; start < kMers.size() && !isStart[start]; ++start);

    std::stringstream superstring;
    superstring << NumberToKMer(kMers[start], k);
    std::vector<bool> mask (1, 1);

    // Move from the first k-mer to the last which has no successor.
    while(edgeFrom[start].secondIndex != size_t(-1)) {
        superstring << NumberToKMer(kMers[edgeFrom[start].secondIndex], k - edgeFrom[start].overlapLength);
        for (int i = 0; i < k - 1 - edgeFrom[start].overlapLength; ++i) mask.push_back(0);
        mask.push_back(1);
        start = edgeFrom[start].secondIndex;
    }

    for(int i = 0; i < k - 1; ++i) mask.emplace_back(0);

    return KMerSet {
            superstring.str(),
            mask,
            k
    };
}

/// Get the approximated shortest superstring of the given k-mers using the GREEDY algorithm.
/// This runs in O(n k), where n is the number of k-mers.
/// If complements are provided, treat k-mer and its complement as identical.
/// If this is the case, k-mers are expected not to contain both k-mer and its complement.
KMerSet Greedy(std::vector<int64_t> kMers, int k, bool complements) {
    if (kMers.empty()) {
        throw std::invalid_argument("input cannot be empty");
    }
    if (complements) {
        size_t n = kMers.size();
        kMers.resize(2 * n);
        for (size_t i = 0; i < n; ++i) {
            kMers[i + n] = ReverseComplement(kMers[i], k);
        }
    }
    auto hamiltonianPath = OverlapHamiltonianPath(kMers, k, complements);
    return SuperstringFromPath(hamiltonianPath, kMers, k);
}
