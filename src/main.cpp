#include "greedy_ac.h"
#include "greedy.h"
#include "generalized_simplitigs.h"
#include "generalized_simplitigs_ac.h"
#include "parser.h"
#include "streaming.h"
#include "output.h"

#include <iostream>
#include <string>
#include "unistd.h"
#include "version.h"

void Help() {
    std::cerr << "KmerCamel version " << VERSION << std::endl;
    std::cerr << "Accepted arguments:" << std::endl;
    std::cerr << "  -p path_to_fasta - required; valid path to fasta file" << std::endl;
    std::cerr << "  -k k_value       - required; integer value for k" << std::endl;
    std::cerr << "  -a algorithm     - the algorithm to be run [global (default), globalAC, local, localAC, streaming]" << std::endl;
    std::cerr << "  -o output_path   - if not specified, the output is printed to stdout" << std::endl;
    std::cerr << "  -d d_value       - integer value for d_max; default 5" << std::endl;
    std::cerr << "  -c               - treat k-mer and its reverse complement as equal" << std::endl;
    std::cerr << "  -h               - print help" << std::endl;
    std::cerr << "  -v               - print version" << std::endl;
    std::cerr << "Example usage:       ./kmercamel -p path_to_fasta -k 13 -d 5 -a local" << std::endl;
    std::cerr << "Possible algorithms: global globalAC local localAC streaming" << std::endl;
}

void Version() {
    std::cerr << VERSION << std::endl;
}

int main(int argc, char **argv) {
    std::string path;
    int k = 0;
    int d_max = 5;
    std::ofstream output;
    std::ostream *of = &std::cout;
    std::string algorithm = "global";
    bool complements = false;
    bool d_set = false;
    int opt;
    try {
        while ((opt = getopt(argc, argv, "p:k:d:a:o:hcv"))  != -1) {
            switch(opt) {
                case  'p':
                    path = optarg;
                    break;
                case 'o':
                    output.open(optarg);
                    of = &output;
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
                case 'v':
                    Version();
                    return 0;
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
        WriteName(k, *of);
        Streaming(path, *of,  k , complements);
    }
    // Handle hash table based separately so that it consumes less memory.
    else if (algorithm == "global" || algorithm == "local") {
        auto kMers = ReadKMers(path, k, complements);
        if (kMers.empty()) {
            std::cerr << "Path '" << path << "' contains no k-mers." << std::endl;
            Help();
            return 1;
        }
        d_max = std::min(k - 1, d_max);
        WriteName(k, *of);
        if (algorithm == "global") Greedy(kMers, *of, k, complements);
        else  GreedyGeneralizedSimplitigs(kMers, *of, k, d_max, complements);
    } else {
        auto data = ReadFasta(path);
        if (data.empty()) {
            std::cerr << "Path '" << path << "' not to a fasta file." << std::endl;
            Help();
            return 1;
        }
        d_max = std::min(k - 1, d_max);

        auto kMers = ConstructKMers(data, k, complements);
        WriteName(k, *of);
        if (algorithm == "globalAC") {
            GreedyAC(kMers, *of, complements);
        }
        else if (algorithm == "localAC") {
            GreedyGeneralizedSimplitigsAC(kMers, *of, k, d_max, complements);
        }
        else {
            std::cerr << "Algorithm '" << algorithm << "' not supported." << std::endl;
            Help();
            return 1;
        }
    }
    *of << std::endl;
    return 0;
}
