#include "greedy.cpp"
#include "generalized_simplitigs.cpp"
#include "generalized_simplitigs_ac.cpp"
#include "kmers.cpp"

#include <iostream>
#include <string>
#include "unistd.h"

#include "bioio.hpp"

void WriteSuperstring(KMerSet result, std::string name) {
    std::cout << name << std::endl;
    std::string superstring = "";
    for (int i = 0; i < result.superstring.length(); ++i) {
        superstring += result.mask[i] ? result.superstring[i] : std::tolower(result.superstring[i]);
    }
    std::cout << superstring << std::endl;
}

void WriteStats(KMerSet result, std::vector<KMer> kMers, std::string data, std::string name) {
    std::cout << "name:                       " << name << std::endl;
    std::cout << "superstring length:         " << result.superstring.length() << std::endl;
    std::cout << "k-mers count:               " << kMers.size() << std::endl;
    std::cout << "length of scanned sequence: " << data.length() << std::endl;
    std::cout << "coefficient:                " << result.superstring.length() / (double)kMers.size() << std::endl;
    std::cout << "========================================="  << std::endl;
}

void Help() {
    std::cerr << "Accepted arguments:" << std::endl;
    std::cerr << "  -p path_to_fasta - required; valid path to fasta file" << std::endl;
    std::cerr << "  -a algortihm     - the algorithm to be run" << std::endl;
    std::cerr << "  -d d_value       - integer value for d_max" << std::endl;
    std::cerr << "  -k k_value       - integer value for k" << std::endl;
    std::cerr << "  -s               - if given print statistics instead of superstring" << std::endl;
    std::cerr << "  -h               - print help" << std::endl;
    std::cerr << "Example usage:       ./kmers -p path_to_fasta -k 13 -d 5 -a greedy" << std::endl;
    std::cerr << "Possible algorithms: greedy pseudosimplitigs pseudosimplitigsAC" << std::endl;
}

int main(int argc, char **argv) {
    std::string path;
    int k = 13;
    int d_max = 5;
    std::string algorithm = "greedy";
    bool printStats = false;
    int opt;
    try {
        while ((opt = getopt(argc, argv, "p:k:d:a:sh"))  != -1) {
            switch(opt) {
                case  'p':
                    path = optarg;
                    break;
                case  'k':
                    k = std::stoi(optarg);
                    break;
                case  'd':
                    d_max = std::stoi(optarg);
                    break;
                case  'a':
                    algorithm = optarg;
                    break;
                case  's':
                    printStats = true;
                    break;
                case 'h':
                    Help();
                    return 0;
            }
        }
    } catch (std::exception) {
        Help();
        return 1;
    };
    auto data = bioio::read_fasta(path);
    if (!data.size()) {
        std::cerr << "Path '" << path << "' not to a fasta file." << std::endl;
        Help();
        return 1;
    }
    for (auto record : data) {
        auto kMers = ConstructKMers(record.sequence, k);
        KMerSet result;
        if (algorithm == "greedy")
            result = Greedy(kMers);
        else if (algorithm == "pseudosimplitigs")
            result = GreedyGeneralizedSimplitigs(kMers, k, d_max);
        else if (algorithm == "pseudosimplitigsAC")
            result = GreedyGeneralizedSimplitigsAC(kMers, k, d_max);
        else {
            std::cerr << "Algortihm '" << algorithm << "' not supported." << std::endl;
            Help();
            return 1;
        }
        if (printStats)
            WriteStats(result, kMers, record.sequence, record.name);
        else
            WriteSuperstring(result, record.name);
        return 0;
    }
}
