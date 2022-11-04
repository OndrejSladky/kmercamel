#include <string>
#include <unordered_set>

#include "models.h"
#include "bioio.hpp"

/// Create a list of unique k-mers in no particular order.
/// This runs in O(k*data.size) expected time.
void AddKMersFromSequence(std::unordered_set<std::string> &kMers, std::string &data, int k) {
    if (size_t(k) > data.size()) return;
    for (size_t i = 0; i <= data.size() - k; ++i) {
        kMers.insert(data.substr(i, k));
    }
}

/// Create a list of unique k-mers in no particular order.
/// This runs in O(k*data.size) expected time.
std::vector<KMer> ConstructKMers(std::vector<bioio::FastaRecord<std::string, std::string>> &data, int k) {
    std::unordered_set<std::string> uniqueKMers;
    for (auto &&record : data) {
        AddKMersFromSequence(uniqueKMers, record.sequence, k);
    }
    std::vector<KMer> result;
    for (auto kMer : uniqueKMers) {
        result.push_back(KMer{kMer});
    }
    return result;
}
