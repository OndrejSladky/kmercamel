#pragma once
#include <string>
#include <vector>
#include <iostream>
#include <queue>
#include <unordered_set>
#include <list>
#include <sstream>
#include <fstream>

#include "models.h"
#include "kmers.h"
#include "global.h"

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
        newState.id = (int)states.size();
        states.push_back(newState);
        return (int)states.size() - 1;
    }

    /// Generate the trie from the given k-mers and set *endStateIndices*.
    void ConstructTrie(const std::vector<KMer> &kMers) {
        endStateIndices = std::vector<int>(kMers.size());
        // Initialize the root.
        AddState(0);
        for (size_t i = 0; i < kMers.size(); ++i) states[0].supporters.push_back(i);
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
        size_t orderingIndex = reversedOrdering.size() - size_t(2);
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

/// Represents oriented edge in the overlap graph.
struct OverlapEdge {
    // Index of the first k-mer.
    size_t firstIndex;
    // Index of the second k-mer.
    size_t secondIndex;
    // Length of the overlap of the two k-mers.
    int overlapLength;
};


/// Greedily find the approximate overlapPath path with longest overlaps using the AC automaton.
std::vector<OverlapEdge> OverlapHamiltonianPathAC (const std::vector<KMer> &kMers, bool complements) {
    size_t n = kMers.size() / (1 + complements);
    ACAutomaton automaton;
    automaton.Construct(kMers);
    std::vector<bool> forbidden(kMers.size(), false);
    std::vector<bool> prefixForbidden(kMers.size(), false);
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
            // If the path forms a cycle, or is between k-mer and its reverse complement, skip this path.
            while (i != incidentKMers[s].end() && (first[*i]%n == j%n || first[*i]%n == last[j]%n|| prefixForbidden[*i])) {
                auto new_i = i;
                new_i++;
                if (prefixForbidden[*i]) incidentKMers[s].erase(i);
                i = new_i;
            }
            if (i == incidentKMers[s].end()) {
                continue;
            }
            std::vector<std::pair<size_t,size_t>> new_edges ({{*i, j}});
            if (complements) new_edges.emplace_back((j + n) % kMers.size(), (*i + n) % kMers.size());
            for (auto [x, y] : new_edges) {
                hamiltonianPath.push_back(OverlapEdge{x, y, automaton.states[s].depth});
                forbidden[y] = true;
                first[last[y]] = first[x];
                last[first[x]] = last[y];
                prefixForbidden[x] = true;
            }
            incidentKMers[s].erase(i);
        }
        incidentKMers[automaton.states[s].backwardEdge].splice(incidentKMers[automaton.states[s].backwardEdge].end(), incidentKMers[s]);
    }
    return hamiltonianPath;
}


/// Return the suffix of the given kMer without the first *overlap* chars.
std::string Suffix(const KMer &kMer, const int overlap) {
    return kMer.value.substr(overlap, kMer.length() - overlap);
}


/// Construct the superstring and the path from the given overlapPath path in the overlap graph.
void SuperstringFromPath(const std::vector<OverlapEdge> &hamiltonianPath, const std::vector<KMer> &kMers, std::ostream& of, const int k) {
    std::vector<OverlapEdge> edgeFrom (kMers.size(), OverlapEdge{size_t(-1),size_t(-1), -1});
    std::vector<bool> isStart(kMers.size(), false);
    // Compute the resulting length of the superstring.
    int result_length = kMers.size() * k;
    for (auto edge : hamiltonianPath) {
        isStart[edge.firstIndex] = true;
        edgeFrom[edge.firstIndex] = edge;
        result_length -= edge.overlapLength;
    }
    for (auto edge : hamiltonianPath) {
        isStart[edge.secondIndex] = false;
    }

    // Find the vertex in the overlap graph with in-degree 0.
    size_t start = 0;
    for (; start < kMers.size() && !isStart[start]; ++start);
    // Handle the edge case of only one k-mer.
    start %= kMers.size();

    // Print the first character.
    of << kMers[start].value[0];

    // Move from the first k-mer to the last which has no successor.
    while(edgeFrom[start].secondIndex != size_t(-1)) {
        int overlapLength = edgeFrom[start].overlapLength;
        if (overlapLength != k - 1) {
            std::string unmaskedNucleotides = kMers[start].value.substr(1, k - 1 - overlapLength);
            std::transform(unmaskedNucleotides.begin(), unmaskedNucleotides.end(), unmaskedNucleotides.begin(), tolower);
            of << unmaskedNucleotides;
        }
        of << kMers[edgeFrom[start].secondIndex].value[0];
        start = edgeFrom[start].secondIndex;
    }

    // Print the trailing k-1 characters.
    std::string unmaskedNucleotides = kMers[start].value.substr(1);
    std::transform(unmaskedNucleotides.begin(), unmaskedNucleotides.end(), unmaskedNucleotides.begin(), tolower);
    of << unmaskedNucleotides;
}

/// Get the approximated shortest superstring of the given k-mers using the global greedy algorithm with Aho-Corasick automaton.
/// This runs in O(n k), where n is the number of k-mers.
/// If complements are provided, it is expected that kMers do not contain both k-mer and its reverse complement.
void GlobalAC(std::vector<KMer> kMers, std::ostream& of, bool complements) {
	if (kMers.empty()) {
		throw std::invalid_argument("input cannot be empty");
	}
    // Add complementary k-mers.
    size_t n = kMers.size();
    kMers.resize(n * (1 + complements));
    if (complements) for (size_t i = 0; i < n; ++i) {
        kMers[i + n] = ReverseComplement(kMers[i]);
    }

	const int k = (int)kMers[0].length();
    auto hamiltonianPath = OverlapHamiltonianPathAC(kMers, complements);
    SuperstringFromPath(hamiltonianPath, kMers, of, k);
}
