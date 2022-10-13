#include <string>
#include <vector>
#include <iostream>
#include <queue>
#include <unordered_set>
#include <list>
#include <sstream>

struct KMer {
	std::string value;

	size_t length() { return value.size(); }
};

struct KMerSet {
	std::string superstring;
	std::vector<bool> mask;
	int k;
};

struct UnionEvent {
	int firstIndex;
	int secondIndex;
	int overlapLength;
};

int NucleotideToInt (char c) {
    switch (c) {
        case 'A': return 0;
        case 'C': return 1;
        case 'G': return 2;
        case 'T': return 3;
        default: throw std::invalid_argument("cannot convert letter " + std::string(1, c) + "to int.");
    }
}

constexpr int INVALID_STATE = -1;
struct ACState {
    std::vector<int> supporters;
    int forwardEdges[4] =  {INVALID_STATE, INVALID_STATE, INVALID_STATE, INVALID_STATE};
    int backwardEdge = 0;
    int depth = 0;
    int id = 0;
};

struct ACAutomaton {
    std::vector<ACState> states;
    std::vector<int> endStateIndices;
    std::vector<int> reversedOrdering;

    int AddState(int depth) {
        ACState newState;
        newState.depth = depth;
        newState.id = states.size();
        states.push_back(newState);
        return states.size() - 1;
    }

    void ConstructTrie(std::vector<KMer> kMers) {
        endStateIndices = std::vector<int>(kMers.size());
        // Initialize the root.
        AddState(0);
        for (int i = 0; i < kMers.size(); ++i) states[0].supporters.push_back(i);
        // TODO: make this independent of the implementation by providing iterator to KMer.
        for (int i = 0; i < kMers.size(); ++i) {
            int state = 0;
            for (char c : kMers[i].value) {
                int index = NucleotideToInt(c);
                if (states[state].forwardEdges[index] == INVALID_STATE) {
                    int newState =  AddState(states[state].depth + 1);
                    states[state].forwardEdges[index] = newState;
                }
                state = states[state].forwardEdges[index];
                states[state].supporters.push_back(i);
            }
            endStateIndices[i] = state;
        }

        for (int i = 0; i < 4; ++i) {
            if (states[0].forwardEdges[i] == INVALID_STATE) {
                states[0].forwardEdges[i] = 0;
            }
        }
    }

    int Step (int state, int index) {
        while (states[state].forwardEdges[index] == INVALID_STATE) {
            state = states[state].backwardEdge;
        }
        return states[state].forwardEdges[index];
    }

    void ConstructBackwardEdges () {
        reversedOrdering = std::vector<int> (states.size(), 0);
        std::queue<int> q;
        for (int i = 0; i < 4; ++i) {
            if (states[0].forwardEdges[i] != 0) {
                q.push(states[0].forwardEdges[i]);
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

    void Construct (std::vector<KMer> kMers) {
        ConstructTrie(kMers);
        ConstructBackwardEdges();
    }
};

struct OverlapEdge {
    int firstIndex;
    int secondIndex;
    int overlapLength;
};

std::vector<OverlapEdge> OverlapHamiltonianPath (std::vector<KMer> &kMers) {
    ACAutomaton automaton;
    automaton.Construct(kMers);
    std::vector<bool> forbidden(kMers.size(), false);
    std::vector<std::list<int>> incidentKMers(kMers.size());
    std::vector<OverlapEdge> hamiltonianPath;
    for (int s : automaton.reversedOrdering) {
        // TODO: implement the kMer algorithm.
    }
    return hamiltonianPath;
}

std::string Suffix(KMer kMer, int overlap) {
    return kMer.value.substr(overlap, kMer.length() - overlap);
}

KMerSet SuperstringFromPath(std::vector<OverlapEdge> &hamiltonianPath, std::vector<KMer> &kMers, int k) {
    std::vector<OverlapEdge> edgeFrom (kMers.size(), OverlapEdge{-1,-1});
    std::vector<bool> isStart(kMers.size(), true);
    for (auto edge : hamiltonianPath) {
        edgeFrom[edge.firstIndex] = edge;
        isStart[edge.secondIndex] = false;
    }

    int start = 0;
    for (; isStart[start] == false; ++start);

    std::stringstream superstring;
    superstring << kMers[start].value;
    std::vector<bool> mask (1, 1);

    while(edgeFrom[start].secondIndex != -1) {
        superstring << Suffix(kMers[edgeFrom[start].secondIndex], edgeFrom[start].overlapLength);
        for (int i = 0; i < k - 1 - edgeFrom[start].overlapLength; ++i) mask.push_back(0);
        mask.push_back(1);
        start = edgeFrom[start].secondIndex;
    }

    for(int i = 0; i < k - 1; ++i) mask.push_back(0);

    return KMerSet {
        superstring.str(),
        mask,
        k
    };
}

/// Get the approximated shortest superstring of the given k-mers using the GREEDY algorithm.
/// This runs in O(n).
KMerSet Greedy(std::vector<KMer> &kMers) {
	if (kMers.size() == 0) {
		throw new std::invalid_argument("input cannot be empty");
	}
	int k = kMers[0].length();
    auto hamiltonianPath = OverlapHamiltonianPath(kMers);
    return SuperstringFromPath(hamiltonianPath, kMers, k);
}
