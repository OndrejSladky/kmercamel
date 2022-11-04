#include <string>
#include <unordered_set>

#include "models.h"

/// Create a list of unique k-mers in no particular order.
/// This runs in O(k*data.size) expected time.
std::vector<KMer> ConstructKMers(std::string data, int k) {
    std::unordered_set<std::string> uniqueKMers;
    if (size_t(k) > data.size()) return {};
    for (size_t i = 0; i <= data.size() - k; ++i) {
        uniqueKMers.insert(data.substr(i, k));
    }
    std::vector<KMer> result;
    for (auto kMer : uniqueKMers) {
        result.push_back(KMer{kMer});
    }
    return result;
}