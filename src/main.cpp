#include <iostream>
#include <string>

#include "unistd.h"
#include "version.h"
#include "ac/global_ac.h"
#include "global.h"
#include "global_sparse.h"
#include "local.h"
#include "ac/local_ac.h"
#include "parser.h"
#include "ac/parser_ac.h"
#include "streaming.h"
#include "khash_utils.h"
#include "conversions.h"

#include <iostream>
#include <string>
#include "unistd.h"
#include "version.h"
#include "masks.h"
#include "lower_bound.h"

int usage() {
    std::cerr << std::endl;
    std::cerr << "Program: KmerCamel - a tool for computation of masked superstrings." << std::endl;
    std::cerr << "Version: " << VERSION << std::endl;
    std::cerr << "Contact: Ondrej Sladky (ondra.sladky@gmail.com)" << std::endl << std::endl;
    std::cerr << "Usage:   kmercamel <command> [options]" << std::endl << std::endl;
    std::cerr << "Command:" << std::endl;
    std::cerr << "    compute    - Compute an approximately shortest masked superstring from input FASTA." << std::endl;
    std::cerr << "    maskopt    - Optimize a given masked superstring." << std::endl;
    std::cerr << "    ms2mssep   - Split a masked superstring into a superstring and a mask." << std::endl;
    std::cerr << "    mssep2ms   - Join a masked superstring from a superstring and a mask." << std::endl;
    std::cerr << "    ms2spss    - Compute rSPSS from a masked superstring (may be of much larger total length)." << std::endl;
    std::cerr << "    spss2ms    - Compute masked superstring corresponding to (r)SPSS." << std::endl;
    std::cerr << "    lowerbound - Compute the lower bound on masked superstring size of a k-mer set." << std::endl;
    std::cerr << std::endl;
    return 1;
}

int usage_subcommand(std::string subcommand) {
    std::cerr << std::endl;
    std::cerr << "Usage:   kmercamel " << subcommand << " [options]";
    if (subcommand == "compute" || subcommand == "lowerbound" || subcommand == "spss2ms")
    std::cerr << " <fasta>";
    else if (subcommand != "mssep2ms")
    std::cerr << " <ms>";

    std::cerr  << std::endl << std::endl;
    std::cerr << "Options:" << std::endl;
    if (subcommand != "mssep2ms" && subcommand != "ms2mssep")
    std::cerr << "  -k INT   - k-mer size [required; up to 127]" << std::endl;

    if (subcommand == "compute")
    std::cerr << "  -a STR   - the algorithm to be run [global (default), streaming, local, globalAC (experimental), localAC (experimental)]" << std::endl;

    else if (subcommand == "maskopt")
    std::cerr << "  -t STR   - the target mask type to be run [maxone (default), minone, minrun, approxminrun]" << std::endl;

    if (subcommand != "lowerbound" && subcommand != "ms2mssep")
    std::cerr << "  -o FILE  - output for the (minone) masked superstring; if not specified, printed to stdout" << std::endl;
    
    if (subcommand == "compute")
    std::cerr << "  -M FILE  - if given, print also ms with mask maximizing ones (only with global)" << std::endl;

    
    if (subcommand == "compute")
    std::cerr << "  -S       - optimize for the input being correctly computed simplitigs (only with global)" << std::endl;

    if (subcommand == "compute")
    std::cerr << "  -d INT   - d_max for local algorithm; default 5" << std::endl;

    if (subcommand == "compute" || subcommand == "maskopt" || subcommand == "lowerbound")
    std::cerr << "  -u       - treat k-mer and its reverse complement as distinct" << std::endl;

    if (subcommand == "mssep2ms") {
    std::cerr << "  -m FILE  - input file with mask" << std::endl;
    std::cerr << "  -s FILE  - input file with superstring" << std::endl;
    }
    if (subcommand == "ms2mssep") {
    std::cerr << "  -m FILE  - output file with mask" << std::endl;
    std::cerr << "  -s FILE  - output file with superstring" << std::endl;
    }
    std::cerr << "  -h       - print help" << std::endl;
    std::cerr << std::endl;
    return 1;
}

constexpr int MAX_K = 127;
constexpr int SIMPLITIG_RATIO_THRESHOLD = 5;

void Version() {
    std::cerr << VERSION << std::endl;
}

/// Run KmerCamel with the given parameters.
template <typename kmer_t, typename kh_wrapper_t>
int kmercamel(kh_wrapper_t wrapper, kmer_t kmer_type, std::string path, int k, int d_max, std::ostream *of, std::ostream *maskf, bool complements, bool masks,
                    std::string algorithm, bool lower_bound, bool assume_simplitigs) {
    if (masks) {
        WriteLog("Started optimization of a masked superstring from '" + path + "'.");
        int ret = Optimize(wrapper, kmer_type, algorithm, path, *of, k, complements);
        if (ret) usage_subcommand("maskopt");
        WriteLog("Finished optimization.");
        return ret;
    }

    if (!lower_bound) WriteLog("Started computation of a masked superstring from '" + path + "'.");
    else WriteLog("Started computation of a masked superstring length lower bound from '" + path + "'.");

    /* Handle streaming algorithm separately. */
    if (algorithm == "streaming") {
        WriteName(path, algorithm, k, false, !complements, *of);
        Streaming(wrapper, kmer_type, path, *of,  k , complements);
        WriteLog("Finished masked superstring computation.");
    }
    /* Handle hash table based algorithms separately so that they consume less memory. */
    else if (algorithm == "global" || algorithm == "local") {

        auto *kMers = wrapper.kh_init_set();
        size_t kmer_count;
        if (!assume_simplitigs) {
            ReadKMers(kMers, wrapper, kmer_type, path, k, complements);

            if (!kh_size(kMers)) {
                std::cerr << "Path '" << path << "' contains no k-mers. Make sure that your file is a FASTA or gzipped FASTA." << std::endl;
                return usage_subcommand("compute");
            }
            kmer_count = kh_size(kMers);
            WriteLog("Finished collecting k-mers: " + std::to_string(kmer_count) + " " + std::to_string(k) + "-mers.");
        }
        
        d_max = std::min(k - 1, d_max);
        if (!lower_bound) WriteName(path, algorithm, k, false, !complements, *of);
        if (maskf != nullptr) WriteName(path, algorithm, k, true, !complements, *maskf);
        if (algorithm == "global") {
            std::vector<simplitig_t> simplitigs;
            if (!assume_simplitigs) {
                simplitigs = get_simplitigs(kMers, wrapper, kmer_type, k, complements);
                wrapper.kh_destroy_set(kMers);
            } else {
                simplitigs = simplitigs_from_fasta(path);
            }
            WriteLog("Finished 1. part: simplitigs (" + std::to_string(simplitigs.size()) + " simplitigs).");
            if (!lower_bound && !assume_simplitigs && simplitigs.size() * SIMPLITIG_RATIO_THRESHOLD >= kmer_count) {
               auto kMerVec = simplitigs_to_kmer_vec(kmer_type, simplitigs, k, kmer_count);
               PartialPreSort(kMerVec, k);
               GlobalSparse(wrapper, kMerVec, *of, maskf, k, complements);
            }
            else if (lower_bound) std::cout << LowerBoundLength(wrapper, kmer_type, simplitigs, k, complements);
            else Global(wrapper, kmer_type, simplitigs, *of, maskf, k, complements);
        } else {
            Local(kMers, wrapper, kmer_type, *of, k, d_max, complements);
            WriteLog("Finished masked superstring computation.");
        }
    } else {
        auto data = ReadFasta(path);
        if (data.empty()) {
            std::cerr << "Path '" << path << "' not to a fasta file." << std::endl;
            return usage_subcommand("compute");
        }
        d_max = std::min(k - 1, d_max);

        auto kMers = ConstructKMers(data, k, complements);
        WriteName(path, algorithm, k, false, !complements, *of);
        if (algorithm == "globalAC") {
            GlobalAC(kMers, *of, complements);
        }
        else if (algorithm == "localAC") {
            LocalAC(kMers, *of, k, d_max, complements);
        }
        else {
            std::cerr << "Algorithm '" << algorithm << "' not supported." << std::endl;
            return usage_subcommand("compute");
        }
        WriteLog("Finished masked superstring computation.");
    }
    *of << std::endl;
    return 0;
}

int camel_compute(int argc, char **argv) {
    std::string subcommand = "compute";
    std::string path;
    if (argc > 1 && std::string(argv[argc - 1]) != "-h") {
        path = argv[argc - 1];
        argc--;
    }
    int k = 0;
    int d_max = 5;
    std::ofstream output;
    std::ostream *of = &std::cout;
    std::ofstream maskOutput;
    std::ostream *maskf = nullptr;
    std::string algorithm = "global";
    bool complements = true;
    bool d_set = false;
    bool assume_simplitigs = false;
    int opt;
    try {
        while ((opt = getopt(argc, argv, "k:d:a:o:huxM:S"))  != -1) {
            switch(opt) {
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
                    break;
                case  'u':
                    complements = false;
                    break;
                case 'x':
                    std::cerr << "Warning: The parameter -x currently has no effect due to the improvement in the underlying algorithm." <<std::endl;
                    break;
                case 'M':
                    maskOutput.open(optarg);
                    maskf = &maskOutput;
                    break;
                case 'S':
                    assume_simplitigs = true;
                    break;
                case 'h':
                    usage_subcommand(subcommand);
                    return 0;
                default:
                    return usage_subcommand(subcommand);
            }
        }
    } catch (std::invalid_argument&) {
        return usage_subcommand(subcommand);
    }
    if (path.empty()) {
        std::cerr << "Required positional parameter path to the file not set." << std::endl;
        return usage_subcommand(subcommand);
    }
    if (k == 0) {
        std::cerr << "Required parameter k not set." << std::endl;
        return usage_subcommand(subcommand);
    } else if (k < 0) {
        std::cerr << "k must be positive." << std::endl;
        return usage_subcommand(subcommand);
    } else if (d_max < 0) {
        std::cerr << "d must be non-negative." << std::endl;
        return usage_subcommand(subcommand);
    } else if (k > MAX_K && (algorithm == "local" || algorithm == "global")) {
        std::cerr << "k > " << MAX_K << " not supported for the algorithm '" + algorithm + "'. Use the  AC version of the algorithm instead." << std::endl;
        return usage_subcommand(subcommand);
    } else if (d_set && (algorithm == "globalAC" || algorithm == "global" || algorithm == "streaming")) {
        std::cerr << "Unsupported argument d for algorithm '" + algorithm + "'." << std::endl;
        return usage_subcommand(subcommand);
    } else if (maskf != nullptr && algorithm != "global") {
        std::cerr << "Outputting mask maximizing number of ones is possible only with global. For other algorithms, use separately kmercamel optimize." << std::endl;
        return usage_subcommand(subcommand);
    } else if (assume_simplitigs && algorithm != "global") {
        std::cerr << "Optimization for the input being simplitigs is possible only with global." << std::endl;
        return usage_subcommand(subcommand);
    }
    if (k < 32) {
        return kmercamel(kmer_dict64_t(), kmer64_t(0), path, k, d_max, of, maskf, complements, false, algorithm, false, assume_simplitigs);
    } else if (k < 64) {
        return kmercamel(kmer_dict128_t(), kmer128_t(0), path, k, d_max, of, maskf, complements, false, algorithm, false, assume_simplitigs);
    } else {
        return kmercamel(kmer_dict256_t(), kmer256_t(0), path, k, d_max, of, maskf, complements, false, algorithm, false, assume_simplitigs);
    }
}

int camel_optimize(int argc, char **argv) {
    std::string subcommand = "maskopt";
    std::string path;
    if (argc > 1 && std::string(argv[argc - 1]) != "-h") {
        path = argv[argc - 1];
        argc--;
    }
    int k = 0;
    std::ofstream output;
    std::ostream *of = &std::cout;
    std::string algorithm = "maxone";
    bool complements = true;
    int opt;
    try {
        while ((opt = getopt(argc, argv, "k:t:o:hu"))  != -1) {
            switch(opt) {
                case 'o':
                    output.open(optarg);
                    of = &output;
                    break;
                case  'k':
                    k = std::stoi(optarg);
                    break;
                case  't':
                    algorithm = optarg;
                    break;
                case  'u':
                    complements = false;
                    break;
                case 'h':
                    usage_subcommand(subcommand);
                    return 0;
                default:
                    return usage_subcommand(subcommand);
            }
        }
    } catch (std::invalid_argument&) {
        return usage_subcommand(subcommand);
    }
    if (path.empty()) {
        std::cerr << "Required positional parameter path to the file not set." << std::endl;
        return usage_subcommand(subcommand);
    }
    if (k == 0) {
        std::cerr << "Required parameter k not set." << std::endl;
        return usage_subcommand(subcommand);
    } else if (k < 0) {
        std::cerr << "k must be positive." << std::endl;
        return usage_subcommand(subcommand);
    } else if (k > MAX_K && (algorithm == "local" || algorithm == "global")) {
        std::cerr << "k > " << MAX_K << " not supported for the algorithm '" + algorithm + "'. Use the  AC version of the algorithm instead." << std::endl;
        return usage_subcommand(subcommand);
    }
    if (k < 32) {
        return kmercamel(kmer_dict64_t(), kmer64_t(0), path, k, 0, of, nullptr, complements, true, algorithm, false, false);
    } else if (k < 64) {
        return kmercamel(kmer_dict128_t(), kmer128_t(0), path, k, 0, of, nullptr, complements, true, algorithm, false, false);
    } else {
        return kmercamel(kmer_dict256_t(), kmer256_t(0), path, k, 0, of, nullptr, complements, true, algorithm, false, false);
    }
}

int camel_lowerbound(int argc, char **argv) {
    std::string subcommand = "lowerbound";
    std::string path;
    if (argc > 1 && std::string(argv[argc - 1]) != "-h") {
        path = argv[argc - 1];
        argc--;
    }
    int k = 0;
    std::ostream *of = &std::cout;
    bool complements = true;
    int opt;
    try {
        while ((opt = getopt(argc, argv, "k:hux"))  != -1) {
            switch(opt) {
                case  'k':
                    k = std::stoi(optarg);
                    break;
                case 'u':
                    complements = false;
                case 'h':
                    usage_subcommand(subcommand);
                    return 0;
                case 'x':
                    std::cerr << "Warning: The parameter -x currently has no effect due to the improvement in the underlying algorithm." <<std::endl;
                    break;
                default:
                    return usage_subcommand(subcommand);
            }
        }
    } catch (std::invalid_argument&) {
        return usage_subcommand(subcommand);
    }
    if (path.empty()) {
        std::cerr << "Required positional parameter path to the file not set." << std::endl;
        return usage_subcommand(subcommand);
    }
    if (k == 0) {
        std::cerr << "Required parameter k not set." << std::endl;
        return usage_subcommand(subcommand);
    } else if (k < 0) {
        std::cerr << "k must be positive." << std::endl;
        return usage_subcommand(subcommand);
    }
    if (k < 32) {
        return kmercamel(kmer_dict64_t(), kmer64_t(0), path, k, 0, of, nullptr, complements, false, "global", true, false);
    } else if (k < 64) {
        return kmercamel(kmer_dict128_t(), kmer128_t(0), path, k, 0, of, nullptr, complements, false, "global", true, false);
    } else {
        return kmercamel(kmer_dict256_t(), kmer256_t(0), path, k, 0, of, nullptr, complements, false, "global", true, false);
    }
}

int camel_split_ms(int argc, char **argv) {
    std::string subcommand = "ms2mssep";
    std::string path;
    if (argc > 1 && std::string(argv[argc - 1]) != "-h") {
        path = argv[argc - 1];
        argc--;
    }

    bool mask_set = false;
    std::ofstream mask;
    std::ostream *maskf = &std::cout;

    bool superstring_set = false;
    std::ofstream superstring;
    std::ostream *superstringf = &std::cout;
    
    int opt;
    try {
        while ((opt = getopt(argc, argv, "m:s:h"))  != -1) {
            switch(opt) {
                case  'm':
                    if (mask_set) {
                        std::cerr << "Error: -m parameter provided multiple times" << std::endl;
                        return usage_subcommand(subcommand);
                    }
                    mask_set = true;
                    mask.open(optarg);
                    maskf = &mask;
                    break;
                case 's':
                    if (superstring_set) {
                        std::cerr << "Error: -s parameter provided multiple times" << std::endl;
                        return usage_subcommand(subcommand);
                    }
                    superstring.open(optarg);
                    superstringf = &superstring;
                    break;
                case 'h':
                    usage_subcommand(subcommand);
                    return 0;
                default:
                    return usage_subcommand(subcommand);
            }
        }
    } catch (std::invalid_argument&) {
        return usage_subcommand(subcommand);
    }
    if (path.empty()) {
        std::cerr << "Required positional parameter path to the file not set." << std::endl;
        return usage_subcommand(subcommand);
    }
    if (!mask_set && !superstring_set) {
        std::cerr << "Cannot have both superstring and mask redirected to stdout." << std::endl;
        return usage_subcommand(subcommand);
    }
    WriteLog("Started splitting masked superstring '" + path + "'.");
    split_ms(*superstringf, *maskf, path);
    WriteLog("Finished masked superstring splitting.");
    return 0;    
}

int camel_join_ms(int argc, char **argv) {
    std::string subcommand = "mssep2ms";

    bool mask_set = false;
    std::ifstream mask;
    std::istream *maskf = &std::cin;

    bool superstring_set = false;
    std::ifstream superstring;
    std::istream *superstringf = &std::cin;

    std::ofstream output;
    std::ostream *of = &std::cout;
    
    
    int opt;
    try {
        while ((opt = getopt(argc, argv, "m:s:o:h"))  != -1) {
            switch(opt) {
                case  'm':
                    if (mask_set) {
                        std::cerr << "Error: -m parameter provided multiple times" << std::endl;
                        return usage_subcommand(subcommand);
                    }
                    mask_set = true;
                    mask.open(optarg);
                    maskf = &mask;
                    break;
                case 's':
                    if (superstring_set) {
                        std::cerr << "Error: -s parameter provided multiple times" << std::endl;
                        return usage_subcommand(subcommand);
                    }
                    superstring_set = true;
                    superstring.open(optarg);
                    superstringf = &superstring;
                    break;
                case 'o':
                    output.open(optarg);
                    of = &output;
                    break;
                case 'h':
                    usage_subcommand(subcommand);
                    return 0;
                default:
                    return usage_subcommand(subcommand);
            }
        }
    } catch (std::invalid_argument&) {
        return usage_subcommand(subcommand);
    }
    if (!mask_set && !superstring_set) {
        std::cerr << "Cannot have both superstring and mask redirected from stdin." << std::endl;
        return usage_subcommand(subcommand);
    }
    WriteLog("Started masked superstring joining.");
    join_ms(*superstringf, *maskf, *of);
    WriteLog("Finished masked superstring joining.");
    return 0;    
}

int camel_ms_to_spss(int argc, char **argv) {
    std::string subcommand = "ms2spss";
    std::string path;
    if (argc > 1 && std::string(argv[argc - 1]) != "-h") {
        path = argv[argc - 1];
        argc--;
    }

    std::ofstream output;
    std::ostream *of = &std::cout;

    int k = 0;
    
    int opt;
    try {
        while ((opt = getopt(argc, argv, "o:k:h"))  != -1) {
            switch(opt) {
                case  'o':
                    output.open(optarg);
                    of = &output;
                    break;
                case 'k':
                    k = std::stoi(optarg);
                    break;
                case 'h':
                    usage_subcommand(subcommand);
                    return 0;
                default:
                    return usage_subcommand(subcommand);
            }
        }
    } catch (std::invalid_argument&) {
        return usage_subcommand(subcommand);
    }
    if (path.empty()) {
        std::cerr << "Required positional parameter path to the file not set." << std::endl;
        return usage_subcommand(subcommand);
    }
    if (k == 0) {
        std::cerr << "Required parameter k not set." << std::endl;
        return usage_subcommand(subcommand);
    } else if (k < 0) {
        std::cerr << "k must be positive." << std::endl;
        return usage_subcommand(subcommand);
    }
    WriteLog("Started rSPSS computation from masked supertring '" + path + "'.");
    ms_to_spss(path, *of, k);
    WriteLog("Finished computing a rSPSS representing the same set.");
    return 0;
}

int camel_spss_to_ms(int argc, char **argv) {
    std::string subcommand = "spss2ms";
    std::string path;
    if (argc > 1 && std::string(argv[argc - 1]) != "-h") {
        path = argv[argc - 1];
        argc--;
    }

    std::ofstream output;
    std::ostream *of = &std::cout;

    int k = 0;
    
    int opt;
    try {
        while ((opt = getopt(argc, argv, "o:k:h"))  != -1) {
            switch(opt) {
                case  'o':
                    output.open(optarg);
                    of = &output;
                    break;
                case 'k':
                    k = std::stoi(optarg);
                    break;
                case 'h':
                    usage_subcommand(subcommand);
                    return 0;
                default:
                    return usage_subcommand(subcommand);
            }
        }
    } catch (std::invalid_argument&) {
        return usage_subcommand(subcommand);
    }
    if (path.empty()) {
        std::cerr << "Required positional parameter path to the file not set." << std::endl;
        return usage_subcommand(subcommand);
    }
    if (k == 0) {
        std::cerr << "Required parameter k not set." << std::endl;
        return usage_subcommand(subcommand);
    } else if (k < 0) {
        std::cerr << "k must be positive." << std::endl;
        return usage_subcommand(subcommand);
    }
    WriteLog("Started masked superstring computation corresponding to (r)SPSS '" + path + "'.");
    spss_to_ms(path, *of, k);
    WriteLog("Finished computing a masked superstring corresponding to the (r)SPSS.");
    return 0;
}

int main(int argc, char **argv) {
    if (argc == 1 || std::string(argv[1]) == "-h") {
        usage();
        return 0;
    } else if (std::string(argv[1]) == "compute") {
        return camel_compute(argc - 1, argv + 1);
    } else if (std::string(argv[1]) == "maskopt") {
        return camel_optimize(argc - 1, argv + 1);
    } else if (std::string(argv[1]) == "lowerbound") {
        return camel_lowerbound(argc - 1, argv + 1);
    } else if (std::string(argv[1]) == "mssep2ms") {
        return camel_join_ms(argc - 1, argv + 1);
    } else if (std::string(argv[1]) == "ms2mssep") {
        return camel_split_ms(argc - 1, argv + 1);
    } else if (std::string(argv[1]) == "ms2spss") {
        return camel_ms_to_spss(argc - 1, argv + 1);
    } else if (std::string(argv[1]) == "spss2ms") {
        return camel_spss_to_ms(argc - 1, argv + 1);
    } else {
        std::cerr << "Command '" << argv[1] << "' not recognized." << std::endl;
        return usage(); 
    }
}
