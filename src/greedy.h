#pragma once
#include "models.h"

#include <vector>
#include <iostream>
#include <unordered_map>
#include <list>
#include <algorithm>
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
    std::vector<size_t> next(kMers.size(), -1);
    for (size_t i = 0; i < kMers.size(); ++i) {
        first[i] = last[i] = i;
    }
    std::unordered_map<int64_t, size_t> prefixes;
    for (int d = k - 1; d >= 0; --d) {
        prefixes.clear();
        for (size_t i = 0 ; i < kMers.size(); ++i) if(!prefixForbidden[i]) {
            next[i] = -1;
            int64_t prefix = BitPrefix(kMers[i], k, d);
            if (prefixes.count(prefix) != 0) next[i] = prefixes[prefix];
            prefixes[prefix] = i;
        }
        for (size_t i = 0 ; i < kMers.size(); ++i) if(!suffixForbidden[i]) {
            int64_t suffix = BitSuffix(kMers[i], d);
            if (prefixes.count(suffix) == 0) continue;
            size_t j = prefixes[suffix];
            size_t previous = prefixes[suffix];
            // If the path forms a cycle, or is between k-mer and its reverse complement, or the k-mers complement was already selected skip this path.
            while (j != size_t(-1) && (first[i]%n == j%n || first[i]%n == last[j]%n || prefixForbidden[j])) {
                size_t new_j = next[j];
                // If the k-mer is forbidden, remove it to keep the complexity linear.
                // This is not done with the first k-mer but that is not a problem.
                if (prefixForbidden[j]) next[previous] = new_j;
                else previous = j;
                j = new_j;
            }
            if (j == size_t(-1)) {
                continue;
            }
            std::vector<std::pair<size_t,size_t>> new_edges ({{i, j}});
            // Add also the edge between complementary k-mers in the opposite direction.
            if (complements) new_edges.emplace_back((j + n) % kMers.size(), (i + n) % kMers.size());
            for (auto [x, y] : new_edges) {
                hamiltonianPath.push_back(OverlapEdge{x, y, d});
                prefixForbidden[y] = true;
                first[last[y]] = first[x];
                last[first[x]] = last[y];
                suffixForbidden[x] = true;
            }
            next[previous] = next[j];
        }
    }
    return hamiltonianPath;
}

/// Construct the superstring and its mask from the given hamiltonian path in the overlap graph.
/// If reverse complements are considered and the hamiltonian path contains two paths which are reverse complements of one another,
/// return only one of them.
void SuperstringFromPath(const std::vector<OverlapEdge> &hamiltonianPath, const std::vector<int64_t> &kMers, std::ostream& of, const int k) {
    std::vector<OverlapEdge> edgeFrom (kMers.size(), OverlapEdge{size_t(-1),size_t(-1), -1});
    std::vector<bool> isStart(kMers.size(), false);
    for (auto edge : hamiltonianPath) {
        isStart[edge.firstIndex] = true;
        edgeFrom[edge.firstIndex] = edge;
    }
    for (auto edge : hamiltonianPath) {
        isStart[edge.secondIndex] = false;
    }

    // Find the vertex in the overlap graph with in-degree 0.
    size_t start = 0;
    for (; start < kMers.size() && !isStart[start]; ++start);

    // Handle the edge case with no edges.
    start %= kMers.size();

    int64_t last = BitSuffix(kMers[start], k-1);
    of << letters[BitPrefix(kMers[start], k, 1)];

    // Move from the first k-mer to the last which has no successor.
    while(edgeFrom[start].secondIndex != size_t(-1)) {
        int overlapLength = edgeFrom[start].overlapLength;
        if (overlapLength != k - 1) {
            std::string unmaskedNucleotides = NumberToKMer(BitPrefix(last, k-1, k-1-overlapLength), k-1-overlapLength);
            std::transform(unmaskedNucleotides.begin(), unmaskedNucleotides.end(), unmaskedNucleotides.begin(), tolower);
            of << unmaskedNucleotides;
        }
        last = BitSuffix(kMers[edgeFrom[start].secondIndex], k-1);
        of << letters[BitPrefix(kMers[edgeFrom[start].secondIndex], k, 1)];
        start = edgeFrom[start].secondIndex;
    }

    // Print the trailing k-1 characters.
    std::string unmaskedNucleotides = NumberToKMer(last, k-1);
    std::transform(unmaskedNucleotides.begin(), unmaskedNucleotides.end(), unmaskedNucleotides.begin(), tolower);
    of << unmaskedNucleotides << "\n";
}

/// Get the approximated shortest superstring of the given k-mers using the GREEDY algorithm.
/// This runs in O(n k), where n is the number of k-mers.
/// If complements are provided, treat k-mer and its complement as identical.
/// If this is the case, k-mers are expected not to contain both k-mer and its complement.
/// Warning: this will destroy kMers.
void Greedy(std::vector<int64_t> &kMers, std::ostream& of, int k, bool complements) {
    if (kMers.empty()) {
        throw std::invalid_argument("input cannot be empty");
    }
    of << ">superstring\n";
    if (complements) {
        size_t n = kMers.size();
        kMers.resize(2 * n);
        for (size_t i = 0; i < n; ++i) {
            kMers[i + n] = ReverseComplement(kMers[i], k);
        }
    }
    auto hamiltonianPath = OverlapHamiltonianPath(kMers, k, complements);
    SuperstringFromPath(hamiltonianPath, kMers, of, k);
}
