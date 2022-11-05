#include <string>
#include <unordered_set>

#include "models.h"
#include "bioio.hpp"

/// Create a list of unique k-mers in no particular order.
/// This runs in O(k*data.size) expected time.
void AddKMersFromSequence(std::unordered_set<std::string> &kMers, std::string data, int k) {
    // Convert the sequence to uppercase letters.
    std::transform(data.begin(), data.end(), data.begin(), toupper);
    for (size_t i = k; i <= data.size(); ++i) {
        if (data[i-1] != 'A' && data[i-1] != 'C' && data[i-1] != 'G' && data[i-1] != 'T') {
            // Skip this and the next k-1 k-mers.
            i += k - 1;
        } else {
            kMers.insert(data.substr(i - k, k));
        }
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
