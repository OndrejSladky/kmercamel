#include "greedy_ac.h"
#include "greedy.h"
#include "generalized_simplitigs.h"
#include "generalized_simplitigs_ac.h"
#include "parser.h"
#include "streaming.h"

#include <iostream>
#include <string>
#include <chrono>
#include "unistd.h"

/// Print the fasta file header.
void WriteName(const KMerSet& result) {
    std::cout << ">superstring ";
    std::cout << "l=" << result.superstring.length() << " ";
    std::cout << "k=" << result.k << std::endl;
}

/// Mask the superstring and print it to the standard output.
void WriteSuperstring(KMerSet result) {
    std::string superstring;
    superstring.reserve(result.superstring.length());
    for (size_t i = 0; i < result.superstring.length(); ++i) {
        superstring +=  result.mask[i] ? (char)result.superstring[i] : (char)std::tolower(result.superstring[i]);
    }
    std::cout << superstring << std::endl;
}

/// Print basic statistics for the resulting superstring.
void WriteStats(const KMerSet& result, long time, size_t kMersCount) {
    std::cout << "superstring length:         " << result.superstring.length() << std::endl;
    std::cout << "number of k-mers:           " << kMersCount << std::endl;
    std::cout << "execution time:             " << time << " ms" << std::endl;
}

void Help() {
    std::cerr << "KmerCamel v0.1" << std::endl;
    std::cerr << "Accepted arguments:" << std::endl;
    std::cerr << "  -p path_to_fasta - required; valid path to fasta file" << std::endl;
    std::cerr << "  -k k_value       - required; integer value for k" << std::endl;
    std::cerr << "  -a algorithm     - the algorithm to be run [global (default), globalAC, local, localAC, streaming]" << std::endl;
    std::cerr << "  -d d_value       - integer value for d_max; default 5" << std::endl;
    std::cerr << "  -s               - if given print statistics instead of superstring" << std::endl;
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
    bool printStats = false;
    bool complements = false;
    bool d_set = false;
    int opt;
    try {
        while ((opt = getopt(argc, argv, "p:k:d:a:shc"))  != -1) {
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
                case  's':
                    printStats = true;
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
    } else if (d_max <= 0) {
        std::cerr << "d must be positive." << std::endl;
        Help();
        return 1;
    } else if (k > 31 && (algorithm == "local" || algorithm == "global")) {
        std::cerr << "k > 31 not supported for the algorithm '" + algorithm + "'. Use its AC version instead." << std::endl;
        Help();
        return 1;
    } else if (d_set && (algorithm == "globalAC" || algorithm == "global")) {
        std::cerr << "Unsupported arguement d for algorithm '" + algorithm + "'." << std::endl;
        Help();
        return 1;
    }

    // Handle streaming algorithm separately and stop the computation afterwards.
    if (algorithm == "streaming") {
        Streaming(path, std::cout,  k , complements);
        return 0;
    }
    auto before = std::chrono::high_resolution_clock::now();
    size_t kMersCount;
    KMerSet result;
    // Handle global separately so that it consumes less memory.
    if (algorithm == "global") {
        auto kMers = ReadKMers(path, k, complements);
        if (kMers.empty()) {
            std::cerr << "Path '" << path << "' contains no k-mers." << std::endl;
            Help();
            return 1;
        }
        kMersCount = kMers.size();
        result = Greedy(kMers, k, complements);
    } else {
        auto data = ReadFasta(path);
        if (data.empty()) {
            std::cerr << "Path '" << path << "' not to a fasta file." << std::endl;
            Help();
            return 1;
        }
        d_max = std::min(k, d_max);

        auto kMers = ConstructKMers(data, k, complements);
        kMersCount = kMers.size();
        if (algorithm == "globalAC")
            result = GreedyAC(kMers, complements);
        else if (algorithm == "local")
            result = GreedyGeneralizedSimplitigs(kMers, k, d_max, complements);
        else if (algorithm == "localAC")
            result = GreedyGeneralizedSimplitigsAC(kMers, k, d_max, complements);
        else {
            std::cerr << "Algorithm '" << algorithm << "' not supported." << std::endl;
            Help();
            return 1;
        }
    }

    auto now = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - before);
    if (printStats) {
        WriteStats(result, duration.count(), kMersCount);
    } else {
        WriteName(result);
        WriteSuperstring(result);
    }
    return 0;
}
