#pragma once
#include <string>
#include <unordered_set>
#include <fstream>
#include <algorithm>

#include <zlib.h>
#include <stdio.h>
#include "kseq.h"
KSEQ_INIT(gzFile, gzread)

#include "kmers.h"
#include "khash_utils.h"
#include <chrono>
#include <iomanip>
#include <ctime>   


/// Fill the k-mer dictionary with k-mers from the given sequence.
/// If complements is true, always add the canonical k-mers.
/// If case_sensitive is true, add the k-mer only if it starts with an upper case letter.
template <typename kmer_t, typename kh_S_t, typename kh_wrapper_t>
void AddKMers(kh_S_t *kMers, kh_wrapper_t wrapper, [[maybe_unused]] kmer_t _, size_t sequence_length,
              const char* sequence, int64_t k, bool complements, bool case_sensitive = false) {
    int64_t currentLength = 0;
    kmer_t currentKMer = 0, reverseComplement = 0;
    kmer_t cases = 0;
    kmer_t mask = (((kmer_t) 1) <<  (2 * k) ) - 1;
    kmer_t shift = 2 * (k - 1);
    for (size_t i = 0; i < sequence_length; ++i) {
        auto data = nucleotideToInt[(uint8_t)sequence[i]];
        if (data >= 4) {
            // Restart if "N"-like nucleotide.
            currentKMer = reverseComplement = 0;
            currentLength = 0;
            continue;
        }
        currentKMer = ((currentKMer << 2) | data) & mask;
        reverseComplement = (reverseComplement >> 2) | ((kmer_t(3 ^ data)) << shift);
        // K-mer is present if it is upper case or case-insensitive.
        cases = (cases | (!case_sensitive || sequence[i] <= 'Z')) << 1;
        if ((++currentLength >= k) && (cases & (kmer_t(1) << k))) {
            // Add the canonical k-mer to the dictionary.
            kmer_t canonical = ((!complements) || currentKMer < reverseComplement) ? currentKMer : reverseComplement;
            int ret;
            wrapper.kh_put_to_set(kMers, canonical, &ret);
        }
    }
}

/// Fill the k-mer dictionary with k-mers from the given sequence and their frequencies.
/// If complements is true, always add the canonical k-mers.
template <typename kmer_t, typename kh_Q_t, typename kh_wrapper_t>
void AddKMersWithFrequencies(kh_Q_t *kMers, kh_wrapper_t wrapper, [[maybe_unused]] kmer_t _, size_t sequence_length,
              const char* sequence, int64_t k, bool complements) {
    int64_t currentLength = 0;
    kmer_t currentKMer = 0, reverseComplement = 0;
    kmer_t mask = (((kmer_t) 1) <<  (2 * k) ) - 1;
    kmer_t shift = 2 * (k - 1);
    for (size_t i = 0; i < sequence_length; ++i) {
        auto data = nucleotideToInt[(uint8_t)sequence[i]];
        if (data >= 4) {
            // Restart if "N"-like nucleotide.
            currentKMer = reverseComplement = 0;
            currentLength = 0;
            continue;
        }
        currentKMer = ((currentKMer << 2) | data) & mask;
        reverseComplement = (reverseComplement >> 2) | ((kmer_t(3 ^ data)) << shift);
        if (++currentLength >= k) {
            // Add the canonical k-mer to the dictionary.
            kmer_t canonical = ((!complements) || currentKMer < reverseComplement) ? currentKMer : reverseComplement;
            auto ptr = wrapper.kh_get_from_freq_map(kMers, canonical);
            bool contained = ptr != kh_end(kMers);
            if (contained) {
                uint16_t count = kh_val(kMers, ptr);
                kh_value(kMers, ptr) = (uint8_t) std::min(255, count + 1);
            } else {
                int ret;
                ptr = wrapper.kh_put_to_freq_map(kMers, canonical, &ret);
                kh_value(kMers, ptr) = 0;
            }
        }
    }
}

/// Return a file/stdin for reading.
gzFile OpenFile(std::string &path) {
    FILE *in_stream;
    if(path=="-"){
        in_stream = stdin;
    }
    else {
        in_stream = fopen(path.c_str(), "r");
        if (in_stream == nullptr) {
            throw std::invalid_argument("couldn't open file " + path);
        }
    }
    gzFile fp = gzdopen(fileno(in_stream), "r");
    return fp;
}


/// Load a dictionary of k-mers from a fasta file.
/// If complements is true, add the canonical k-mers.
template <typename kmer_t, typename kh_S_t, typename kh_wrapper_t>
void ReadKMers(kh_S_t *kMers, kh_wrapper_t wrapper, kmer_t _, std::string &path, int k, bool complements,
               bool case_sensitive = false) {
    gzFile fp = OpenFile(path);
    kseq_t *seq = kseq_init(fp);

    while (kseq_read(seq) >= 0) {
        AddKMers(kMers, wrapper, _, seq->seq.l, seq->seq.s, k, complements, case_sensitive);
    }

    kseq_destroy(seq);
    gzclose(fp);
}

/// Load a dictionary of k-mers from a fasta file.
/// If complements is true, add the canonical k-mers.
template <typename kmer_t, typename kh_S_t, typename kh_wrapper_t>
void ReadKMersFiltered(kh_S_t *kMers, kh_wrapper_t wrapper, kmer_t _, std::string &path, int k, bool complements,
               uint16_t min_frequency) {
    gzFile fp = OpenFile(path);
    kseq_t *seq = kseq_init(fp);
    
    auto kMersFrequencies = wrapper.kh_init_freq_map();

    while (kseq_read(seq) >= 0) {
        AddKMersWithFrequencies(kMersFrequencies, wrapper, _, seq->seq.l, seq->seq.s, k, complements);
    }


    auto kMersVec = kMersToVecFiltered(kMersFrequencies, _, min_frequency);
    wrapper.kh_destroy_freq_map(kMersFrequencies);
    kMersFromVec(kMers, wrapper, kMersVec);

    kseq_destroy(seq);
    gzclose(fp);
}


/// Read the masked superstring from the given path and return it wrapped as a kseq_t.
kseq_t* ReadMaskedSuperstring(std::string &path) {
    gzFile fp = OpenFile(path);
    kseq_t *seq = kseq_init(fp);
    kseq_read(seq);
    return seq;
}

/// Ensure that the file is at the end.
void AssertEOF(kseq_t *seq, std::string message) {
    if (kseq_read(seq) >= 0) {
        throw std::invalid_argument(message);
    }
}

void WriteLog(const std::string message) {
    auto snapshot = std::chrono::system_clock::now();
    std::time_t time = std::chrono::system_clock::to_time_t(snapshot);
    std::tm* time_tm = std::localtime(&time);
    std::cerr << "[" << std::put_time(time_tm, "%H:%M:%S") << "] " << message << std::endl;
}

/// Print the fasta file header.
void WriteName(const std::string &dataset, const std::string &algorithm, const int k, const bool maxone, const bool unidirectional, std::ostream &of) {
    of << ">maskedsuperstring ";
    of << "dataset='" << dataset << "' ";
    of << "k=" << k << " ";
    of << "alg=" << algorithm << " ";
    of << "mask=";
    if (maxone) of << "maxone ";
    else of << "minone ";
    of << "mode=";
    if (unidirectional) of << "unidirectional";
    else of << "bidirectional";
    of << std::endl;
}
