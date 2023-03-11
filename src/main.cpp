#include "greedy_ac.h"
#include "greedy.h"
#include "generalized_simplitigs.h"
#include "generalized_simplitigs_ac.h"
#include "parser.h"
#include "streaming.h"
#include "output.h"

#include <iostream>
#include <string>
#include <chrono>
#include "unistd.h"

void Help() {
    std::cerr << "KmerCamel v0.1" << std::endl;
    std::cerr << "Accepted arguments:" << std::endl;
    std::cerr << "  -p path_to_fasta - required; valid path to fasta file" << std::endl;
    std::cerr << "  -k k_value       - required; integer value for k" << std::endl;
    std::cerr << "  -a algorithm     - the algorithm to be run [global (default), globalAC, local, localAC, streaming]" << std::endl;
    std::cerr << "  -d d_value       - integer value for d_max; default 5" << std::endl;
    std::cerr << "  -c               - treat k-mer and its reverse complement as equal" << std::endl;
    std::cerr << "  -h               - print help" << std::endl;
    std::cerr << "Example usage:       ./kmers -p path_to_fasta -k 13 -d 5 -a global" << std::endl;
    std::cerr << "Possible algorithms: global globalAC local localAC streaming" << std::endl;
}

int main(int argc, char **argv) {
    std::string path;
    int k = 0;
    int d_max = 5;
    std::string algorithm = "global";
    bool complements = false;
    bool d_set = false;
    int opt;
    try {
        while ((opt = getopt(argc, argv, "p:k:d:a:hc"))  != -1) {
            switch(opt) {
                case  'p':
                    path = optarg;
                    break;
                case  'k':
                    k = std::stoi(optarg);
                    break;
                case  'd':
                    d_set = true;
                    d_max = std::stoi(optarg);
                    break;
                case  'a':
                    algorithm = optarg;
                    // Backwards compatability.
                    if (algorithm == "greedy") algorithm = "global";
                    if (algorithm == "greedyAC") algorithm = "globalAC";
                    if (algorithm == "pseudosimplitigs") algorithm = "local";
                    if (algorithm == "pseudosimplitigsAC") algorithm = "localAC";
                    break;
                case  'c':
                    complements = true;
                    break;
                case 'h':
                default:
                    Help();
                    return 0;
            }
        }
    } catch (std::invalid_argument&) {
        Help();
        return 1;
    }
    if (path.empty()) {
        std::cerr << "Required parameter p not set." << std::endl;
        Help();
        return 1;
    }
    if (k == 0) {
        std::cerr << "Required parameter k not set." << std::endl;
        Help();
        return 1;
    } else if (k < 0) {
        std::cerr << "k must be positive." << std::endl;
        Help();
        return 1;
    } else if (d_max < 0) {
        std::cerr << "d must be non-negative." << std::endl;
        Help();
        return 1;
    } else if (k > 31 && (algorithm == "local" || algorithm == "global")) {
        std::cerr << "k > 31 not supported for the algorithm '" + algorithm + "'. Use its AC version instead." << std::endl;
        Help();
        return 1;
    } else if (d_set && (algorithm == "globalAC" || algorithm == "global" || algorithm == "streaming")) {
        std::cerr << "Unsupported arguement d for algorithm '" + algorithm + "'." << std::endl;
        Help();
        return 1;
    }

    // Handle streaming algorithm separately.
    if (algorithm == "streaming") {
        WriteName(k);
        Streaming(path, std::cout,  k , complements);
    }
    // Handle greedy separately so that it consumes less memory.
    else if (algorithm == "global") {
        auto kMers = ReadKMers(path, k, complements);
        if (kMers.empty()) {
            std::cerr << "Path '" << path << "' contains no k-mers." << std::endl;
            Help();
            return 1;
        }
        WriteName(k);
        Greedy(kMers, std::cout, k, complements);
    } else {
        auto data = ReadFasta(path);
        if (data.empty()) {
            std::cerr << "Path '" << path << "' not to a fasta file." << std::endl;
            Help();
            return 1;
        }
        d_max = std::min(k - 1, d_max);

        auto kMers = ConstructKMers(data, k, complements);
        WriteName(k);
        if (algorithm == "globalAC") {
            KMerSet result = GreedyAC(kMers, std::cout, complements);
            WriteSuperstring(result.superstring, result.mask);
        }
        else if (algorithm == "local") {
            GreedyGeneralizedSimplitigs(kMers, std::cout, k, d_max, complements);
        }
        else if (algorithm == "localAC") {
            GreedyGeneralizedSimplitigsAC(kMers, std::cout, k, d_max, complements);
        }
        else {
            std::cerr << "Algorithm '" << algorithm << "' not supported." << std::endl;
            Help();
            return 1;
        }
    }
    std::cout << std::endl;
    return 0;
}
