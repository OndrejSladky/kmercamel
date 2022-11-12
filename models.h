#pragma once
#include <string>
#include <vector>

struct KMer {
    std::string value;
    size_t length(){ return value.size(); }
};

struct KMerSet {
    std::string superstring;
    std::vector<bool> mask;
    int k;
};
