#pragma once
#include "models.h"
#include "greedy_ac.h"

#include <string>
#include <vector>
#include <deque>
#include <list>

/// Find the index of the first extending k-mer from the incidentKMers which is not forbidden.
/// Mark this k-mer forbidden and remove all the k-mers in incidentKMers which are forbidden.
/// The latter is done in order not to increase time complexity by repeatedly iterating over forbidden k-mers.
/// Return -1 if only forbidden k-mers were found.
size_t ExtensionAC(std::vector<bool> &forbidden, std::list<size_t> &incidentKMers) {
    while(!incidentKMers.empty()) {
        if(forbidden[incidentKMers.front()]) {
            incidentKMers.pop_front();
        } else {
            forbidden[incidentKMers.front()] = true;
            size_t ret = incidentKMers.front();
            incidentKMers.pop_front();
            return ret;
        }
    }
    return -1;
}

/// Compute the generalized simplitigs greedily using the Aho-Corasick automaton.
/// This runs in O(n k), where n is the number of k-mers.
KMerSet GreedyGeneralizedSimplitigsAC(std::vector<KMer> kMers, int k, int d_max) {
    std::string superstring = "";
    std::vector<bool> mask;
    ACAutomaton a;
    a.Construct(kMers);

    // suffixes[i][j] is the state in the AC automaton given by the suffix of kMers[i] of size j (or -1 is none).
    std::vector<std::vector<int>> suffixes(kMers.size(), std::vector<int> (k + 1, -1));
    // suffixes[i][j] is the state in the AC automaton given by the prefix of kMers[i] of size j.
    std::vector<std::vector<int>> prefixes(kMers.size(), std::vector<int>(k + 1, 0));
    // For each state the list of k-mers which have the given state as a suffix.
    std::vector<std::list<size_t>> incidentKMers(a.states.size());
    // true if the given k-mer has already been used.
    std::vector<bool> forbidden(kMers.size(), false);

    for (size_t i = 0; i < kMers.size(); ++i) {
        for(int j = 0; j < k; ++j) {
            prefixes[i][j + 1] = a.states[prefixes[i][j]].forwardEdges[NucleotideToInt(kMers[i].value[j])];
        }
        for (int s = a.endStateIndices[i]; ; s = a.states[s].backwardEdge) {
            suffixes[i][a.states[s].depth] = s;
            incidentKMers[s].push_back(i);
            if (s == 0) break;
        }
    }

    size_t firstUnused = 0;
    for(;;) {
        // Find the first unused k-mer.
        while(forbidden[firstUnused]) {
            ++firstUnused;
            if (firstUnused == kMers.size()) {
                firstUnused = size_t(-1);
                break;
            }
        }
        if (firstUnused == size_t(-1)) break;
        auto simplitig = kMers[firstUnused].value;
        // Maintain the left and right most k-mer of the generalized simplitig.
        int firstKMer = firstUnused;
        int lastKMer = firstUnused;
        forbidden[firstUnused] = true;
        std::deque<bool> simplitigMask{1};
        int d_l = 1, d_r = 1;
        while (d_l <= d_max || d_r <= d_max) {
            if (d_r <= d_l) {
                int state = suffixes[lastKMer][k - d_r];
                int ext = -1;
                if (state != -1) ext = ExtensionAC(forbidden, a.states[state].supporters);
                if (ext == -1) {
                    // No right extension found.
                    ++d_r;
                } else {
                    // Extend the simplitig to the right.
                    lastKMer = ext;
                    simplitig += kMers[ext].value.substr(k - d_r, d_r);
                    for (int i = 0; i < d_r - 1; ++i) simplitigMask.push_back(0);
                    simplitigMask.push_back(1);
                    d_r = 1;
                }
            } else {
                int state = prefixes[firstKMer][k - d_l];
                int ext = -1;
                if (state != -1) ext = ExtensionAC(forbidden, incidentKMers[state]);
                if (ext == -1) {
                    // No left extension found.
                    ++d_l;
                } else {
                    // Extend the simplitig to the left.
                    firstKMer = ext;
                    simplitig = kMers[ext].value.substr(0, d_l) + simplitig;
                    for (int i = 0; i < d_l - 1; ++i) simplitigMask.push_front(0);
                    simplitigMask.push_front(1);
                    d_l = 1;
                }
            }
        }
        superstring += simplitig;
        for (auto x: simplitigMask) mask.push_back(x);
        for (int i = 0; i < k - 1; ++i) mask.push_back(0);
    }

    return KMerSet{superstring, mask, k};
}

