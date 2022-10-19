#include "greedy.cpp"
#include "kmers.cpp"

#include <iostream>
#include <string>

#include "bioio.hpp"


void WriteResult(KMerSet result, std::vector<KMer> kMers, std::string data, std::string name) {
    std::cout << "name:                       " << name << std::endl;
    std::cout << "superstring length:         " << result.superstring.length() << std::endl;
    std::cout << "k-mers count:               " << kMers.size() << std::endl;
    std::cout << "length of scanned sequence: " << data.length() << std::endl;
    std::cout << "coefficient:                " << result.superstring.length() / (double)kMers.size() << std::endl;
    std::cout << "========================================="  << std::endl;
}

int main(int argc, char **argv) {
    std::string path;
    int k;
    try {
        path = argv[1];
        k = std::stoi(argv[2]);
    } catch (std::exception) {
        std::cerr << "Expected usage: ./kmers path_to_fasta k" << std::endl;
        return 1;
    };
    auto data = bioio::read_fasta(path);
    for (auto record : data) {
        auto kMers = ConstructKMers(record.sequence, k);
        auto result = Greedy(kMers);
        WriteResult(result, kMers, record.sequence, record.name);
        return 0;
    }
}
