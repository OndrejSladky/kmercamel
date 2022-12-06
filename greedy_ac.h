#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <queue>
#include <unordered_set>
#include <list>
#include <sstream>

#include "models.h"
#include "kmers.h"
#include "greedy.h"

constexpr int INVALID_STATE = -1;
struct ACState {
    // List of k-mers whose prefix this state is.
    std::list<size_t> supporters;
    // Indexes of the states longer by the corresponding nucleotide.
    int forwardEdges[4] =  {INVALID_STATE, INVALID_STATE, INVALID_STATE, INVALID_STATE};
    // Where to go if searching failed.
    int backwardEdge = 0;
    // Length of the corresponding string.
    int depth = 0;
    // Index of the state in the automaton's *states* list.
    int id = 0;
};

struct ACAutomaton {
    std::vector<ACState> states;
    // Indices of the states where i-th k-mer ends.
    std::vector<int> endStateIndices;
    // Reversed BFS order.
    std::vector<int> reversedOrdering;

    /// Append a new state and return its ID.
    int AddState(const int depth) {
        ACState newState;
        newState.depth = depth;
        newState.id = states.size();
        states.push_back(newState);
        return states.size() - 1;
    }

    /// Generate the trie from the given k-mers and set *endStateIndices*.
    void ConstructTrie(const std::vector<KMer> &kMers) {
        endStateIndices = std::vector<int>(kMers.size());
        // Initialize the root.
        AddState(0);
        for (size_t i = 0; i < kMers.size(); ++i) states[0].supporters.push_back(i);
        // TODO: make this independent of the implementation by providing iterator to KMer.
        for (size_t i = 0; i < kMers.size(); ++i) {
            int state = 0;
            for (char c : kMers[i].value) {
                int index = NucleotideToInt(c);
                // If the next state does not yet exist, create it.
                if (states[state].forwardEdges[index] == INVALID_STATE) {
                    int newState =  AddState(states[state].depth + 1);
                    states[state].forwardEdges[index] = newState;
                }
                state = states[state].forwardEdges[index];
                states[state].supporters.push_back(i);
            }
            endStateIndices[i] = state;
        }

        // Create a forward edge from the root to itself so that the AC Step always finds a valid forward edge.
        for (int & forwardEdge : states[0].forwardEdges) {
            if (forwardEdge == INVALID_STATE) {
                forwardEdge = 0;
            }
        }
    }

    /// Do one step of the AC algorithm from the given state with a nucleotide with a given index.
    int Step (int state, const int index) {
        while (states[state].forwardEdges[index] == INVALID_STATE) {
            state = states[state].backwardEdge;
        }
        return states[state].forwardEdges[index];
    }

    /// Construct the fail edges for the trie that already had been created.
    /// Also construct the reversed BFS ordering.
    void ConstructBackwardEdges () {
        reversedOrdering = std::vector<int> (states.size(), 0);
        std::queue<int> q;
        for (int & forwardEdge : states[0].forwardEdges) {
            if (forwardEdge != 0) {
                q.push(forwardEdge);
            }
        }
        int orderingIndex = reversedOrdering.size() - 2;
        while (!q.empty()) {
            int state = q.front();
            reversedOrdering[orderingIndex--] = state;
            q.pop();
            for (int i = 0; i < 4; ++i) {
                int nextState = states[state].forwardEdges[i];
                if (nextState != INVALID_STATE && nextState != 0) {
                    states[nextState].backwardEdge = Step(states[state].backwardEdge, i);
                    q.push(nextState);
                }
            }
        }
    }

    /// Construct the Aho-Corasick automaton.
    void Construct (const std::vector<KMer> &kMers) {
        ConstructTrie(kMers);
        ConstructBackwardEdges();
    }
};

/// Greedily find the approximate hamiltonian path with longest overlaps using the AC automaton.
std::vector<OverlapEdge> OverlapHamiltonianPathAC (const std::vector<KMer> &kMers) {
    ACAutomaton automaton;
    automaton.Construct(kMers);
    std::vector<bool> forbidden(kMers.size(), false);
    std::vector<std::list<size_t>> incidentKMers(automaton.states.size());
    std::vector<OverlapEdge> hamiltonianPath;
    std::vector<size_t> first(kMers.size());
    std::vector<size_t> last(kMers.size());
    for (size_t i = 0; i < kMers.size(); ++i) {
        first[i] = last[i] = i;
        incidentKMers[automaton.states[automaton.endStateIndices[i]].backwardEdge].push_back(i);
    }
    for (int s : automaton.reversedOrdering) {
        if (incidentKMers[s].empty()) continue;
        for (size_t j : automaton.states[s].supporters) {
            if (incidentKMers[s].empty()) continue;
            if (forbidden[j]) continue;
            auto i = incidentKMers[s].begin();
            if (first[*i] == j) {
                if (incidentKMers[s].size() == 1) continue;
                i++;
            }
            hamiltonianPath.push_back(OverlapEdge{*i, j, automaton.states[s].depth});
            forbidden[j] = true;
            first[last[j]] = first[*i];
            last[first[*i]] = last[j];
            incidentKMers[s].erase(i);
        }
        incidentKMers[automaton.states[s].backwardEdge].merge(incidentKMers[s]);
    }
    return hamiltonianPath;
}

/// Return the suffix of the given kMer without the first *overlap* chars.
std::string Suffix(const KMer &kMer, const int overlap) {
    return kMer.value.substr(overlap, kMer.length() - overlap);
}

/// Construct the superstring from the given hamiltonian path in the overlap graph.
KMerSet SuperstringFromPath(const std::vector<OverlapEdge> &hamiltonianPath, const std::vector<KMer> &kMers, const int k) {
    std::vector<int64_t> encoded(kMers.size());
    for (size_t i = 0 ; i < kMers.size(); ++i) {
        encoded[i] = KMerToNumber(kMers[i]);
    }
    return SuperstringFromPath(hamiltonianPath, encoded, k);
}

/// Get the approximated shortest superstring of the given k-mers using the GREEDY algorithm with Aho-Corasick automaton.
/// This runs in O(n k), where n is the number of k-mers.
KMerSet GreedyAC(std::vector<KMer> &kMers) {
	if (kMers.empty()) {
		throw std::invalid_argument("input cannot be empty");
	}
	const int k = kMers[0].length();
    auto hamiltonianPath = OverlapHamiltonianPathAC(kMers);
    return SuperstringFromPath(hamiltonianPath, kMers, k);
}
