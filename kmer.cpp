#include <algorithm>

#include "utils.hpp"
#include "kmer.hpp"

read_kmers_t extract_kmers_from_read(std::string read, int kmer_size, bool both_strands) {
    std::string rev_read = reverse_complement(read);

    std::vector<kmer_t> read_kmers = std::vector<kmer_t>(read.length()-kmer_size);
    std::vector<kmer_t> rev_read_kmers = std::vector<kmer_t>(rev_read.length()-kmer_size);

    kmer_bv_t read_bv_kmers;
    kmer_bv_t rev_read_bv_kmers;

    // CHECK: if its worth having everything in 1 for, or it can be split in 2 functions
    // without losing performance (for the sake of clarity)
    for (int k = 0; k < read.length() - std::min(KMER_BV_SIZE, kmer_size); ++k) {
        if (k < read.length() - kmer_size) {
            auto hash = hash_kmer(read.substr(k, kmer_size));
            read_kmers[k] = kmer_t(hash, k);

            if (both_strands) {
                auto rev_hash = hash_kmer(rev_read.substr(k, kmer_size));
                rev_read_kmers[k] = kmer_t(rev_hash, k);
            }
        }

        if (k < read.length() - KMER_BV_SIZE) {
            auto hash = hash_kmer(read.substr(k, KMER_BV_SIZE));
            read_bv_kmers.set(hash);
            
            if (both_strands) {
                auto rev_hash = hash_kmer(rev_read.substr(k, KMER_BV_SIZE));
                rev_read_bv_kmers.set(rev_hash);
            }
        }
    }

    std::sort(read_kmers.begin(), read_kmers.end());
    if (both_strands) std::sort(rev_read_kmers.begin(), rev_read_kmers.end());
    return read_kmers_t{read_kmers, rev_read_kmers, read_bv_kmers, rev_read_bv_kmers};
}

read_kmers_t extract_minimizers_from_read(std::string read, int kmer_size, int window_size, bool both_strands) {
    std::string rev_read = reverse_complement(read);

    std::vector<kmer_t> read_kmers = std::vector<kmer_t>();
    std::vector<kmer_t> rev_read_kmers = std::vector<kmer_t>();

    kmer_bv_t read_bv_kmers;
    kmer_bv_t rev_read_bv_kmers;

    for (int p = 0; p < read.length() - kmer_size - window_size; ++p) {
        auto all_windows = read.substr(p, kmer_size+window_size-1);
        auto all_windows_rev = rev_read.substr(p, kmer_size+window_size-1);

        int32_t min_kmer = -1;
        int32_t rev_min_kmer = -1;

        int32_t pos = 0;
        int32_t pos_rev = 0;

        for (int w = 0; w < window_size; ++w) {
            auto hash = hash_kmer(all_windows.substr(w, kmer_size));
            if (min_kmer == -1 || hash < min_kmer) {
                min_kmer = hash;
                pos = w;
            }

            if (both_strands) {
                auto rev_hash = hash_kmer(all_windows_rev.substr(w, kmer_size));
                if (rev_min_kmer == -1 || rev_hash < rev_min_kmer) {
                    rev_min_kmer = rev_hash;
                    pos_rev = w;
                }
            }
        }

        if (read_kmers.size() == 0 || read_kmers[read_kmers.size()-1].second != p+pos) read_kmers.push_back(kmer_t(min_kmer, p+pos));
        if (both_strands && (rev_read_kmers.size() == 0 || rev_read_kmers[rev_read_kmers.size()-1].second != p+pos_rev)) rev_read_kmers.push_back(kmer_t(rev_min_kmer, p+pos_rev));
    }

    std::sort(read_kmers.begin(), read_kmers.end());
    if (both_strands) std::sort(rev_read_kmers.begin(), rev_read_kmers.end());
    return read_kmers_t{read_kmers, rev_read_kmers, read_bv_kmers, rev_read_bv_kmers};
}

// k1 and k2 are sorted by pair first (kmer hash)
std::vector<kmer_match_t> get_common_kmers(const std::vector<kmer_t> &k1, const std::vector<kmer_t> &k2) {
    int p1 = k1.size() - 1;
    int p2 = k2.size() - 1;
    std::vector<kmer_match_t> intersection;
    intersection.reserve(p1 + 1);

    while (p1 >= 0 && p2 >= 0) {
        while (p2 >= 0 && k2[p2].first > k1[p1].first) {
            --p2;
        }

        int p2t = p2;
        while (p2 >= 0 && k2[p2].first == k1[p1].first) {
            intersection.push_back(kmer_match_t(k1[p1].second, k2[p2].second));
            --p2;
        }
        p2 = p2t;
        --p1;
    }

    std::sort(intersection.begin(), intersection.end());
    return intersection;
}