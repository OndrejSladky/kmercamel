#include "models.h"
#include "greedy.cpp"

#include <string>
#include <vector>
#include <deque>
#include <list>

int ExtensionAC(std::vector<bool> &forbidden, std::list<int> &incidentKMers) {
    while(!incidentKMers.empty()) {
        if(forbidden[incidentKMers.front()]) {
            incidentKMers.pop_front();
        } else {
            forbidden[incidentKMers.front()] = true;
            int ret = incidentKMers.front();
            incidentKMers.pop_front();
            return ret;
        }
    }
    return -1;
}

KMerSet GreedyGeneralizedSimplitigsAC(std::vector<KMer> kMers, int k, int d_max) {
    std::string superstring = "";
    std::vector<bool> mask;
    ACAutomaton a;
    a.Construct(kMers);

    std::vector<std::vector<int>> suffixes(kMers.size(), std::vector<int> (k + 1, -1));
    std::vector<std::vector<int>> prefixes(kMers.size(), std::vector<int>(k + 1, 0));
    std::vector<std::list<int>> incidentKMers(a.states.size());
    std::vector<bool> forbidden(kMers.size(), false);

    for (int i = 0; i < kMers.size(); ++i) {
        for(int j = 0; j < k; ++j) {
            prefixes[i][j + 1] = a.states[prefixes[i][j]].forwardEdges[NucleotideToInt(kMers[i].value[j])];
        }
        for (int s = a.endStateIndices[i]; ; s = a.states[s].backwardEdge) {
            suffixes[i][a.states[s].depth] = s;
            incidentKMers[s].push_back(i);
            if (s == 0) break;
        }
    }

    int firstUnused = 0;

    for(;;) {
        while(forbidden[firstUnused]) {
            ++firstUnused;
            if (firstUnused == kMers.size()) {
                firstUnused = -1;
                break;
            }
        }
        if (firstUnused == -1) break;
        auto simplitig = kMers[firstUnused].value;
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
                    ++d_r;
                } else {
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
                    ++d_l;
                } else {
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

