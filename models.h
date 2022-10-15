#ifndef KMERS_MODELS_H
#define KMERS_MODELS_H
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

#endif //KMERS_MODELS_H
