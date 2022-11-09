#include "greedy_ac.cpp"
#include "greedy.cpp"
#include "generalized_simplitigs.cpp"
#include "generalized_simplitigs_ac.cpp"
#include "parser.cpp"

#include <iostream>
#include <string>
#include <chrono>
#include <vector>
#include "unistd.h"


void WriteSuperstring(KMerSet result) {
    std::string superstring = "";
    for (size_t i = 0; i < result.superstring.length(); ++i) {
        superstring += result.mask[i] ? result.superstring[i] : std::tolower(result.superstring[i]);
    }
    std::cout << superstring << std::endl;
}

void WriteStats(KMerSet result, std::vector<KMer> kMers, std::vector<FastaRecord> &data, long time) {
    size_t scanned_length = 0;
    for (auto &&record : data) {
        scanned_length += record.sequence.length();
    }
    std::cout << "superstring length:         " << result.superstring.length() << std::endl;
    std::cout << "k-mers count:               " << kMers.size() << std::endl;
    std::cout << "length of scanned sequence: " << scanned_length << std::endl;
    std::cout << "coefficient:                " << result.superstring.length() / (double)kMers.size() << std::endl;
    std::cout << "execution time:             " << time << " ms" << std::endl;
    std::cout << "========================================="  << std::endl;
}

void Help() {
    std::cerr << "Accepted arguments:" << std::endl;
    std::cerr << "  -p path_to_fasta - required; valid path to fasta file" << std::endl;
    std::cerr << "  -k k_value       - required; integer value for k" << std::endl;
    std::cerr << "  -a algortihm     - the algorithm to be run [greedy (default), greedyAC, pseudosimplitigs, pseudosimplitigsAC]" << std::endl;
    std::cerr << "  -d d_value       - integer value for d_max; default 5" << std::endl;
    std::cerr << "  -s               - if given print statistics instead of superstring" << std::endl;
    std::cerr << "  -c               - treat k-mer and its reverse complement as equal" << std::endl;
    std::cerr << "  -h               - print help" << std::endl;
    std::cerr << "Example usage:       ./kmers -p path_to_fasta -k 13 -d 5 -a greedy" << std::endl;
    std::cerr << "Possible algorithms: greedy greedyAC pseudosimplitigs pseudosimplitigsAC" << std::endl;
}

int main(int argc, char **argv) {
    std::string path;
    int k = 0;
    int d_max = 5;
    std::string algorithm = "greedy";
    bool printStats = false;
    bool complements = false;
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
                    d_max = std::stoi(optarg);
                    break;
                case  'a':
                    algorithm = optarg;
                    break;
                case  's':
                    printStats = true;
                    break;
                case  'c':
                    complements = true;
                    break;
                case 'h':
                    Help();
                    return 0;
            }
        }
    } catch (std::invalid_argument) {
        Help();
        return 1;
    };
    auto data = ReadFasta(path);
    if (!data.size()) {
        std::cerr << "Path '" << path << "' not to a fasta file." << std::endl;
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
    } else if (k > 31 && (algorithm == "greedy" || algorithm == "pseudosimplitigs")) {
        std::cerr << "k > 31 not supported for the algorithm '" + algorithm + "'. Use its AC version instead." << std::endl;
        Help();
        return 1;
    }

    auto kMers = ConstructKMers(data, k);
    auto before = std::chrono::high_resolution_clock::now();
    KMerSet result;
    if (algorithm == "greedyAC")
        result = GreedyAC(kMers);
    else if (algorithm == "greedy")
        result = Greedy(kMers);
    else if (algorithm == "pseudosimplitigs")
        result = GreedyGeneralizedSimplitigs(kMers, k, d_max, complements);
    else if (algorithm == "pseudosimplitigsAC")
        result = GreedyGeneralizedSimplitigsAC(kMers, k, d_max);
    else {
        std::cerr << "Algortihm '" << algorithm << "' not supported." << std::endl;
        Help();
        return 1;
    }
    auto now = std::chrono::high_resolution_clock::now();
    if (printStats) {
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - before);
        WriteStats(result, kMers, data, duration.count());
    } else {
        WriteSuperstring(result);
    }
    return 0;
}
