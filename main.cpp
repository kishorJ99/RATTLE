#include "fasta.hpp"
#include "cluster.hpp"
#include "utils.hpp"
#include "correct.hpp"
#include "argagg.hpp"
#include "kmer.hpp"
#include "hps/src/hps.h"
#include "spoa/spoa.hpp"
// #include "DBSCAN/clustering.cpp"
#include "SimpleDBSCAN/dbscan.h"

#include <iostream>
#include <future>
#include <fstream>
#include <stdio.h>
#include <algorithm>
#include <queue>
#include <unistd.h>
#include <cmath>

read_set_t read_multiple_inputs(std::vector<std::string> input_files, std::vector<std::string> label_files, bool raw, int lower_len, int upper_len) {
        
        read_set_t reads;

        bool no_labels = label_files.size() == 0;

        if (input_files.size() != label_files.size() && !no_labels) {
            throw "\nError: Number of input files and number of label files do not match\n";
        }

        int sample_number = 0;

        for (std::string i : input_files) {
            if(access(i.c_str(), F_OK )){
                throw "\nError: Input file not found! \n";
            } else {

                std::string sample_label = no_labels ? "" : "." + label_files[sample_number];
                std::string filename = i;
                int index = filename.find_last_of(".");
                std::string extension = filename.substr(index + 1);

                if (!extension.compare("gz")){
                    filename = unzip_file(filename, index);
                    index = filename.find_last_of(".");
                    extension = filename.substr(index + 1);
                }

                if (!extension.compare("fq") || !extension.compare("fastq")){
                    auto file_reads = read_fastq_file(filename, sample_label, raw, lower_len, upper_len);

                    reads.insert(std::end(reads), std::begin(file_reads), std::end(file_reads));

                } else if (!extension.compare("fasta") || !extension.compare("fa")){
                    auto file_reads = read_fasta_file(filename, sample_label, raw, lower_len, upper_len);
                    reads.insert(std::end(reads), std::begin(file_reads), std::end(file_reads));
                } else {
                    throw "\nError: Input file format incorrect! Please use fasta/fastq file. \n";
                }

                ++sample_number;
            }
        }

        return reads;
}

std::vector<std::string> splitString(std::string str, char delimiter) {
    std::vector<std::string> internal;
    std::stringstream ss(str); // Turn the string into a stream.
    std::string tok;

    while(getline(ss, tok, delimiter)) {
        internal.push_back(tok);
    }

    return internal;
}

struct vec1s {
    std::string data[1];
    std::string sample;
    std::string operator[](int idx) const { return data[idx]; }
};
struct vec1f {
    float data[1];
    float operator[](int idx) const { return data[idx]; }
};

struct vec4096f {
    float data[4096];
    float operator[](int idx) const { return data[idx]; }
};


float getDistanceFunction(vec1s p1, vec1s p2, bool is_rna) {

    // Trim reads
    std::string p1Trimmed;
    std::string p2Trimmed;

    if (p1[0].size() > p2[0].size()) {
        p1Trimmed = p1[0].substr(p1[0].size() - p2[0].size(), p2[0].size());
        p2Trimmed = p2[0];
    } else {
        p1Trimmed = p1[0];
        p2Trimmed = p2[0].substr(p2[0].size() - p1[0].size(), p1[0].size());
    }

    // Choose Optomal Values.
    int kmerSize = ceil((log(p1Trimmed.size()) / log(4))) ;
    float errorRate = 0.05;
    float radius = ceil(errorRate * p1Trimmed.size() * kmerSize) ;



    // Generate Kmers
    kmer_bv k1 = extract_kmers_method_2(kmerSize, p1Trimmed);
    kmer_bv k2 = extract_kmers_method_2(kmerSize, p2Trimmed);
    
    float distance = static_cast<float>((k1 ^ k2).count()) / static_cast<float>((k1 | k2).count());

    // float test = 0.54320;

    std::cerr << distance << " Top : " <<  (k1 ^ k2).count() << " Bottom : " << ((k1 | k2).count()) << std::endl;

    std::ofstream outfile;
    outfile.open("Distance Matrix.csv", std::ios_base::app);

    outfile << p1.sample << ", " << p2.sample <<  ", " << radius << ", " << kmerSize << ", " << distance << ", " << p1Trimmed.size() << ", \n";

    outfile.close();

    if (p1.sample == p2.sample) {
        std::cerr << "K1 : " << p1.sample << " K2 : " << p2.sample <<  " Radius : " << radius << " Kmer Size : " << kmerSize << " Distance : " << distance << " Read Length : " << p1Trimmed.size() << std::endl;
    }
    // std::cerr << "K1 : " << k1BV << std::endl;
    // std::cerr << "K2 : " << k2BV << std::endl;
    return distance;
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cout << "Run with mode: ./rattle <cluster|cluster_summary|extract_clusters|correct|polish>" << std::endl;
        return EXIT_FAILURE;
    }

    char* mode = argv[1];
    if (!strcmp(mode, "cluster")) {
        argagg::parser argparser {{
            { "help", {"-h", "--help"},
            "shows this help message", 0},
            { "input", {"-i", "--input"},
            "input fasta/fastq file (required)", 1},
            { "label", {"-l", "--label"},
            "labels for the files in order of entry", 1},
            { "output", {"-o", "--output"},
            "output folder (default: .)", 1},
            { "threads", {"-t", "--threads"},
            "number of threads to use (default: 1)", 1},
            { "kmer_size", {"-k", "--kmer-size"},
            "k-mer size for gene clustering (default: 10, maximum: 16)", 1},
            { "t_s", {"-s", "--score-threshold"},
            "minimum score for two reads to be in the same gene cluster (default: 0.2)", 1},  
            { "t_v", {"-v", "--max-variance"},
            "max allowed variance for two reads to be in the same gene cluster (default: 1000000)", 1},
            { "iso", {"--iso"},
            "perform clustering at the isoform level", 0},
            { "iso_kmer_size", {"--iso-kmer-size"},
            "k-mer size for isoform clustering (default: 11, maximum: 16)", 1},
            { "iso_t_s", {"--iso-score-threshold"},
            "minimum score for two reads to be in the same isoform cluster (default: 0.3)", 1},  
            { "iso_t_v", {"--iso-max-variance"},
            "max allowed variance for two reads to be in the same isoform cluster (default: 25)", 1},
            { "bv_threshold", {"-B", "--bv-start-threshold"},
            "starting threshold for the bitvector k-mer comparison (default: 0.4)", 1},  
            { "bv_min_threshold", {"-b", "--bv-end-threshold"},
            "ending threshold for the bitvector k-mer comparison (default: 0.2)", 1},  
            { "bv_falloff", {"-f", "--bv-falloff"},
            "falloff value for the bitvector threshold for each iteration (default: 0.05)", 1},  
            { "min_reads_cluster", {"-r", "--min-reads-cluster"},
            "minimum number of reads per cluster (default: 0)", 1},
            { "repr_percentile", {"-p", "--repr-percentile"},
            "cluster representative percentile (default: 0.15)", 1},  
            { "rna", {"--rna"},
            "use this mode if data is direct RNA (disables checking both strands)", 0},
            { "verbose", {"--verbose"},
            "use this flag if need to print the progress", 0},
            { "raw", {"--raw"},
            "use this flag if want to use raw datasets", 0},
            {"lower_len", {"--lower-length"},
            "set the lower length for input reads filter (default: 150)", 1},
            {"upper_len", {"--upper-length"},
            "set the upper length for input reads filter (default: 100,000)", 1},
            {"iso_freq", {"--iso-freq"},
            "use this mode to cluster iso forms by length.", 0},
            {"iso_bv", {"--iso-bitvec"},
            "use this mode to cluster iso forms by bit vectors.", 0},
            {"iso_len_bv", {"--iso-len-bv"},
            "use this mode to cluster iso forms by length then bit vectors.", 0},
            {"dbscan_min_points", {"--dbscan-mp"},
            "Use this to pass in DBscan parameters 'min points'", 1},
            {"dbscan_eps", {"--dbscan-eps"},
            "Use this to pass in DBscan parameter 'Eps'", 1},
            {"iso_trim", {"--iso-trim"},
            "Use this to choose trimming reads before creating bv", 0}
        }};

        argagg::parser_results args;
        try {
            args = argparser.parse(argc, argv);
        } catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
            return EXIT_FAILURE;
        }

        if (args["help"]) {
            std::cerr << argparser;
            return EXIT_SUCCESS;
        }

        if (!args["input"]) {
            std::cerr << "ERROR: No input file provided" << std::endl;
            std::cerr << argparser;
            return EXIT_FAILURE;
        }

        int n_threads = args["threads"].as<int>(1);

        int kmer_size = args["kmer_size"].as<int>(10);
        double t_s = args["t_s"].as<double>(0.2);
        double t_v = args["t_v"].as<double>(1000000);

        int iso_kmer_size = args["iso_kmer_size"].as<int>(11);
        double iso_t_s = args["iso_t_s"].as<double>(0.3);
        double iso_t_v = args["iso_t_v"].as<double>(25);

        double bv_threshold = args["bv_threshold"].as<double>(0.4);
        double bv_min_threshold = args["bv_min_threshold"].as<double>(0.2);
        double bv_falloff = args["bv_falloff"].as<double>(0.05);

        int min_reads_cluster = args["min_reads_cluster"].as<int>(0);
        double repr_percentile = args["repr_percentile"].as<double>(0.15);

        int lower_len = args["lower_len"].as<int>(150);
        int upper_len = args["upper_len"].as<int>(100000);

        bool verbose = args["verbose"];
        bool raw = args["raw"];

        if(kmer_size > 16 || iso_kmer_size > 16){
            std::cerr << "\nError: maximum kmer size = 16 \n";
            return EXIT_FAILURE;
        }

        bool is_rna = args["rna"];

        int db_scan_mp = args["dbscan_min_points"].as<int>(3);
        float db_scan_eps = args["dbscan_eps"].as<float>(5);

        std::cerr << "RNA mode: " << std::boolalpha << is_rna << std::endl;

        std::cerr << "Reading fasta file... " << std::endl;

        read_set_t reads;

        std::vector<std::string> files = splitString(args["input"].as<std::string>(""), ',');
        std::vector<std::string> labels = splitString(args["label"].as<std::string>(""), ',');;

        try {
           reads = read_multiple_inputs(files, labels, raw, lower_len, upper_len);
        }
        catch (char* c) {
            std::cerr << c;
            return EXIT_FAILURE;
        }

        std::cout << "Reads: " << reads.size() << std::endl;

        sort_read_set(reads);

        std::cerr << "Done" << std::endl;

        auto gene_clusters = cluster_reads(reads, kmer_size, t_s, t_v, bv_threshold, bv_min_threshold, bv_falloff, min_reads_cluster, false, repr_percentile, is_rna, verbose, n_threads, false);
        std::ofstream out_file(args["output"].as<std::string>(".") + "/clusters.out", std::ofstream::binary);

        std::cerr << "Gene clustering done" << std::endl;
        std::cerr << gene_clusters.size() << " gene clusters found" << std::endl;
        if (!args["iso"] && !args["iso_freq"] && !args["iso_bv"] && !args["iso_len_bv"]) {
            hps::to_stream(gene_clusters, out_file);
            out_file.close();
            return EXIT_SUCCESS;
        }

        // clustering at isoform level
        cluster_set_t iso_clusters;

        if (args["iso"]) {
            int i = 0;
            for (auto &c : gene_clusters) {
                // sort gene cluster seqs by size
                std::stable_sort(c.seqs.begin(), c.seqs.end(), [&reads](cseq_t a, cseq_t b) {
                    return a.seq_id > b.seq_id;
                });

                std::stable_sort(c.seqs.begin(), c.seqs.end(), [&reads](cseq_t a, cseq_t b) {
                    return reads[a.seq_id].seq.size() > reads[b.seq_id].seq.size();
                });

                // generate new read set with gene cluster reads
                read_set_t gene_reads;
                for (auto &cs : c.seqs) {
                    gene_reads.push_back(reads[cs.seq_id]);
                }

                bool iso_trim = args['iso_trim'];

                // cluster gene reads & save new iso clusters
                auto iso_clusters_tmp = cluster_reads(gene_reads, iso_kmer_size, iso_t_s, iso_t_v, bv_threshold, bv_min_threshold, bv_falloff, min_reads_cluster, false, repr_percentile, is_rna, false, n_threads, iso_trim);
                for (auto &ic : iso_clusters_tmp) {
                    cluster_t iso_cluster;
                    iso_cluster.main_seq = cseq_t{c.seqs[ic.main_seq.seq_id].seq_id, ic.main_seq.rev};

                    for (auto &ics : ic.seqs) {
                        iso_cluster.seqs.push_back(cseq_t{c.seqs[ics.seq_id].seq_id, ics.rev});
                    }

                    iso_clusters.push_back(iso_cluster);
                }

                ++i;
                if (verbose) print_progress(i, gene_clusters.size());
            }
        } else if (args["iso_freq"]) {

            int distributionThreshold = args["iso_freq"].as<int>(100);

            int i = 0;
            for (auto &c : gene_clusters) {

                // create point vector
                std::vector<vec1f> points;

                for (auto &cs : c.seqs) {
                    points.push_back(vec1f{(float) reads[cs.seq_id].seq.size()});
                }

                auto dbscan = DBSCAN<vec1f, float>();

                dbscan.Run(&points, 1, 2.0f, 3);

                auto noise = dbscan.Noise;
                auto clusters = dbscan.Clusters;
                
                // convert clusters to iso clusters

                for (auto &cluster : clusters) {
                    cluster_t iso_cluster;

                    iso_cluster.main_seq = cseq_t{c.seqs[cluster[0]].seq_id, c.seqs[cluster[0]].rev};

                    for (auto &cs : cluster) {
                        iso_cluster.seqs.push_back(cseq_t{c.seqs[cs].seq_id, c.seqs[cs].rev});
                    }

                    iso_clusters.push_back(iso_cluster);
                }
                ++i;
            }
        } else if (args["iso_bv"]) {
            std::cerr << "Mode Bitvec" << std::endl;
            std::cerr << "eps : " << db_scan_eps << " mp : " << db_scan_mp << " k : " << iso_kmer_size << std::endl;

            // extract kmers from reads
            std::vector<std::vector<kmer_t>> kmers(reads.size());
            std::vector<std::vector<kmer_t>> rev_kmers(reads.size());

            std::vector<kmer_bv_t> bv_kmers(reads.size());
            std::vector<kmer_bv_t> rev_bv_kmers(reads.size());
            
            std::vector<std::future<void>> tasks;
            for (int t = 0; t < n_threads; ++t) {
                tasks.emplace_back(std::async(std::launch::async, [t, &reads, n_threads, iso_kmer_size, &kmers, &rev_kmers, &bv_kmers, &rev_bv_kmers, is_rna] {
                    for (int i = t; i < reads.size(); i+=n_threads) {
                        read_kmers_t k1 = extract_kmers_from_read(reads[i].seq, iso_kmer_size, !is_rna);

                        kmers[i] = k1.list_forward;
                        rev_kmers[i] = k1.list_reverse;
                        bv_kmers[i] = k1.bv_forward;
                        rev_bv_kmers[i] = k1.bv_reverse;
                    }
                }));
            }

            for (auto &&task : tasks) {
                task.get();
            }

            int i=0;
            for (auto &c : gene_clusters) {

                // create point vector
                std::vector<vec4096f> points;

                // bv outputfile
                std::ofstream file;
                file.open ("example_" + std::to_string(i) + ".csv");

                for (auto &cs : c.seqs) {
                    vec4096f vect;
                    std::string fileLine = "";
                    auto bv = bv_kmers[cs.seq_id];

                    for (int i = 0; i < 4096; ++i) {
                        vect.data[i] = (float) bv[i];
                        fileLine += std::to_string(bv[i]) + ",";
                    }
                    fileLine += "\n";
 
                    file << fileLine;
                    points.push_back(vect);
                }

                file.close();

                auto dbscan = DBSCAN<vec4096f, float>();

                dbscan.Run(&points, 4096, db_scan_eps, db_scan_mp);

                auto noise = dbscan.Noise;
                auto clusters = dbscan.Clusters;
                
                // std::cerr << "Cluster Size : " <<clusters.size() << std::endl;

                // convert clusters to iso clusters
                int j = 0;
                for (auto &cluster : clusters) {
                    cluster_t iso_cluster;

                    iso_cluster.main_seq = cseq_t{c.seqs[cluster[0]].seq_id, c.seqs[cluster[0]].rev};

                    for (auto &cs : cluster) {
                        iso_cluster.seqs.push_back(cseq_t{c.seqs[cs].seq_id, c.seqs[cs].rev});
                    }

                    iso_clusters.push_back(iso_cluster);
                }
                ++j;
                ++i;
            }
        } else if (args["iso_len_bv"]) {
            std::cerr << "Mode length bitvec" << std::endl;
            std::cerr << "eps : " << db_scan_eps << " mp : " << db_scan_mp << " k : " << iso_kmer_size << std::endl;

            int i=0;
            for (auto &c : gene_clusters) {

                // input vector
                std::vector<vec1s> points; 

                for (auto  &cs : c.seqs) {
                    vec1s vect;

                    vect.data[0] = reads[cs.seq_id].seq;
                    vect.sample = reads[cs.seq_id].header.substr(reads[cs.seq_id].header.find_last_of(".")) + ":" + std::to_string(cs.seq_id);

                    points.push_back(vect);
                }

                std::function<float(vec1s p1, vec1s p2)> distanceFunction = is_rna ? [](vec1s p1, vec1s p2) {return getDistanceFunction(p1, p2, true);} : [](vec1s p1, vec1s p2) {return getDistanceFunction(p1, p2, false);};
                
                auto dbscan = DBSCAN<vec1s, float>();
                dbscan.Run(&points, 1, 10, 2, distanceFunction);

                auto noise = dbscan.Noise;
                auto clusters = dbscan.Clusters;
                
                std::cerr << "No of Clusters : " <<clusters.size() << " Cluster Size : " << clusters[0].size() << std::endl;

                // if (clusters.size() > 1) {
                //     // convert clusters to iso clusters
                //     int j = 0;
                    for (auto &cluster : clusters) {
                        cluster_t iso_cluster;

                        iso_cluster.main_seq = cseq_t{c.seqs[cluster[0]].seq_id, c.seqs[cluster[0]].rev};

                        for (auto &cs : cluster) {
                            iso_cluster.seqs.push_back(cseq_t{c.seqs[cs].seq_id, c.seqs[cs].rev});
                        }

                        iso_clusters.push_back(iso_cluster);
                    }
                //     ++j;
                //     ++i;
                // } else {
                    // iso_clusters.push_back(c);
                // }
            }

        }


        std::cerr << "Isoform clustering done" << std::endl;
        std::cerr << iso_clusters.size() << " isoform clusters found" << std::endl;
        hps::to_stream(iso_clusters, out_file);
        out_file.close();
        return EXIT_SUCCESS;
    } else if (!strcmp(mode, "correct")) {
        argagg::parser argparser {{
            { "help", {"-h", "--help"},
            "shows this help message", 0},
            { "input", {"-i", "--input"},
            "input fasta/fastq file (required)", 1},
            { "label", {"-l", "--label"},
            "labels for the files in order of entry", 1},
            { "clusters", {"-c", "--clusters"},
            "clusters file (required)", 1},
            { "output", {"-o", "--output"},
            "output folder (default: .)", 1},
            { "gap-occ", {"-g", "--gap-occ"},
            "gap-occ (default: 0.3)", 1},
            { "min-occ", {"-m", "--min-occ"},
            "min-occ (default: 0.3)", 1},
            { "split", {"-s", "--split"},
            "split clusters into sub-clusters of size s for msa (default: 200)", 1},
            { "min-reads", {"-r", "--min-reads"},
            "min reads to correct/output consensus for a cluster (default: 5)", 1},
            { "threads", {"-t", "--threads"},
            "number of threads to use (default: 1)", 1},
            { "verbose", {"--verbose"},
            "use this flag if need to print the progress", 0},
            {"raw", {"--raw"},
            "use this flag if want to use raw datasets", 0},
            {"lower_len", {"--lower-length"},
            "set the lower length for input reads filter (default: 150)", 1},
            {"upper_len", {"--upper-length"},
            "set the upper length for input reads filter (default: 100,000)", 1},
        }};

        argagg::parser_results args;
        try {
            args = argparser.parse(argc, argv);
        } catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
            return EXIT_FAILURE;
        }

        if (args["help"]) {
            std::cerr << argparser;
            return EXIT_SUCCESS;
        }

        if (!args["input"]) {
            std::cerr << "ERROR: No input file provided" << std::endl;
            std::cerr << argparser;
            return EXIT_FAILURE;
        }

        if (!args["clusters"]) {
            std::cerr << "ERROR: No clusters file provided" << std::endl;
            std::cerr << argparser;
            return EXIT_FAILURE;
        }

        bool raw = args["raw"];
        int lower_len = args["lower_len"].as<int>(150);
        int upper_len = args["upper_len"].as<int>(100000);

        std::cerr << "Reading fasta file... ";

        read_set_t reads;

        std::vector<std::string> files = splitString(args["input"].as<std::string>(""), ',');
        std::vector<std::string> labels = splitString(args["label"].as<std::string>(""), ',');;

        try {
           reads = read_multiple_inputs(files, labels, raw, lower_len, upper_len);
        }
        catch (char* c) {
            std::cerr << c;
            return EXIT_FAILURE;
        }
           
        sort_read_set(reads);
        std::cerr << "Done" << std::endl;

        int n_threads = args["threads"].as<int>(1);
        std::ifstream in_file(args["clusters"].as<std::string>(), std::ifstream::binary);
        auto clusters = hps::from_stream<cluster_set_t>(in_file);
        int split = args["split"].as<int>(200);
        double min_occ = args["min-occ"].as<double>(0.3);
        double gap_occ = args["gap-occ"].as<double>(0.3);
        int min_reads = args["min-reads"].as<int>(5);
        bool verbose = args["verbose"];

        correction_results_t correction = correct_reads(clusters, reads, min_occ, gap_occ, 30.0, split, min_reads, n_threads, verbose);
        write_fastq_file(correction.corrected, args["output"].as<std::string>(".") + "/corrected.fq");
        write_fastq_file(correction.uncorrected, args["output"].as<std::string>(".") + "/uncorrected.fq");
        write_fastq_file(correction.consensi, args["output"].as<std::string>(".") + "/consensi.fq");

        std::cerr << "Done" << std::endl;
        
        // std::cout << alignmentGraph << std::endl;
    } else if (!strcmp(mode, "cluster_summary")) {
        argagg::parser argparser {{
            { "help", {"-h", "--help"},
            "shows this help message", 0},
            { "input", {"-i", "--input"},
            "input fasta/fastq file (required)", 1},
            { "label", {"-l", "--label"},
            "labels for the files in order of entry", 1},
            { "clusters", {"-c", "--clusters"},
            "clusters file (required)", 1},
            { "raw", {"--raw"},
            "use this flag if want to use raw datasets", 0},  
            { "lower_len", {"--lower-length"},
            "set the lower length for input reads filter (default: 150)", 1},
            { "upper_len", {"--upper-length"},
            "set the upper length for input reads filter (default: 100,000)", 1},         
        }};

        argagg::parser_results args;
        try {
            args = argparser.parse(argc, argv);
        } catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
            return EXIT_FAILURE;
        }

        if (args["help"]) {
            std::cerr << argparser;
            return EXIT_SUCCESS;
        }

        if (!args["input"]) {
            std::cerr << "ERROR: No input file provided" << std::endl;
            std::cerr << argparser;
            return EXIT_FAILURE;
        }

        if (!args["clusters"]) {
            std::cerr << "ERROR: No clusters file provided" << std::endl;
            std::cerr << argparser;
            return EXIT_FAILURE;
        }

        std::cerr << "Reading fasta file... ";
        
        bool raw = args["raw"];
        int lower_len = args["lower_len"].as<int>(150);
        int upper_len = args["upper_len"].as<int>(100000);
        read_set_t reads;

        std::vector<std::string> files = splitString(args["input"].as<std::string>(""), ',');
        std::vector<std::string> labels = splitString(args["label"].as<std::string>(""), ',');;

        try {
           reads = read_multiple_inputs(files, labels, raw, lower_len, upper_len);
        }
        catch (char* c) {
            std::cerr << c;
            return EXIT_FAILURE;
        }

        sort_read_set(reads);
        std::cerr << "Done" << std::endl;

        std::ifstream in_file(args["clusters"].as<std::string>(), std::ifstream::binary);
        auto clusters = hps::from_stream<cluster_set_t>(in_file);

        int cid = 0;
        for (auto c : clusters) {
            for (auto seq : c.seqs) {
                std::cout << reads[seq.seq_id].header << "," << cid << std::endl;
            }

            ++cid;
        }
    } else if (!strcmp(mode, "extract_clusters")) {
        argagg::parser argparser {{
            { "help", {"-h", "--help"},
            "shows this help message", 0},
            { "input", {"-i", "--input"},
            "input fasta/fastq file (required)", 1},
            { "label", {"-l", "--label"},
            "labels for the files in order of entry", 1},
            { "clusters", {"-c", "--clusters"},
            "clusters file (required)", 1},
            { "output", {"-o", "--output-folder"},
            "output folder for fastx files (default: .)", 1},
            { "minreads", {"-m", "--min-reads"},
            "min reads per cluster to save it into a file", 1},
            { "fastq", {"--fastq"},
            "whether input and output should be in fastq format (instead of fasta)", 0},
            {"raw", {"--raw"},
            "use this flag if want to use raw datasets", 0},         
            { "lower_len", {"--lower-length"},
            "set the lower length for input reads filter (default: 150)", 1},
            { "upper_len", {"--upper-length"},
            "set the upper length for input reads filter (default: 100,000)", 1},
        }};

        argagg::parser_results args;
        try {
            args = argparser.parse(argc, argv);
        } catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
            return EXIT_FAILURE;
        }

        if (args["help"]) {
            std::cerr << argparser;
            return EXIT_SUCCESS;
        }

        if (!args["input"]) {
            std::cerr << "ERROR: No input file provided" << std::endl;
            std::cerr << argparser;
            return EXIT_FAILURE;
        }

        if (!args["clusters"]) {
            std::cerr << "ERROR: No clusters file provided" << std::endl;
            std::cerr << argparser;
            return EXIT_FAILURE;
        }

        std::cerr << "Reading fasta file... ";
        
        bool raw = args["raw"];
        int lower_len = args["lower_len"].as<int>(150);
        int upper_len = args["upper_len"].as<int>(100000);
        read_set_t reads;

        std::vector<std::string> files = splitString(args["input"].as<std::string>(""), ',');
        std::vector<std::string> labels = splitString(args["label"].as<std::string>(""), ',');;

        try {
           reads = read_multiple_inputs(files, labels, raw, lower_len, upper_len);
        }
        catch (char* c) {
            std::cerr << c;
            return EXIT_FAILURE;
        }

        sort_read_set(reads);
        std::cerr << "Done" << std::endl;

        std::ifstream in_file(args["clusters"].as<std::string>(), std::ifstream::binary);
        auto clusters = hps::from_stream<cluster_set_t>(in_file);
        int min_reads = args["minreads"].as<int>(0);

        int cid = 0;
        for (auto c : clusters) {
            if (c.seqs.size() > min_reads) {
                std::ostringstream ss_fn;
                if (args["output"]) {
                    ss_fn << args["output"].as<std::string>();
                    ss_fn << "/";
                }

                ss_fn << "cluster_";
                ss_fn << cid;

                if (args["fastq"]) {
                    ss_fn << ".fq";
                } else {
                    ss_fn << ".fa";
                }

                std::ofstream cfile;
                cfile.open(ss_fn.str());
                
                for (auto seq : c.seqs) {
                    // std::cout << reads[seq.seq_id].header << "," << cid << std::endl;
                    cfile << reads[seq.seq_id].header << "\n";
                    if (seq.rev) {
                        cfile << reverse_complement(reads[seq.seq_id].seq) << "\n";
                    } else {
                        cfile << reads[seq.seq_id].seq << "\n";
                    }

                    if (args["fastq"]) {
                        cfile << reads[seq.seq_id].ann << "\n";
                        cfile << reads[seq.seq_id].quality << "\n";
                    }
                }
                
                cfile.close();
            }

            ++cid;
        }
    } else if (!strcmp(mode, "polish")) {
        argagg::parser argparser {{
            { "help", {"-h", "--help"},
            "shows this help message", 0},
            { "input", {"-i", "--input"},
            "input RATTLE consensi fasta/fastq file (required)", 1},
            { "output", {"-o", "--output-folder"},
            "output folder for fastx files (default: .)", 1},
            { "threads", {"-t", "--threads"},
            "number of threads to use (default: 1)", 1},
            { "rna", {"--rna"},
            "use this mode if data is direct RNA (disables checking both strands)", 0},
            { "verbose", {"--verbose"},
            "use this flag if need to print the progress", 0},
        }};

        argagg::parser_results args;
        try {
            args = argparser.parse(argc, argv);
        } catch (const std::exception& e) {
            std::cerr << e.what() << std::endl;
            return EXIT_FAILURE;
        }

        if (args["help"]) {
            std::cerr << argparser;
            return EXIT_SUCCESS;
        }

        if (!args["input"]) {
            std::cerr << "ERROR: No input file provided" << std::endl;
            std::cerr << argparser;
            return EXIT_FAILURE;
        }

        std::cerr << "Reading fasta file... ";
        if(access(args["input"], F_OK )){
            std::cerr << "\nError: Input file not found! \n";
            return EXIT_FAILURE;
        }
        
        read_set_t reads = read_fastq_file(args["input"], true, 150, 100000);

        sort_read_set(reads);
        std::cerr << "Done" << std::endl;

        int n_threads = args["threads"].as<int>(1);
        bool is_rna = args["rna"];
        bool verbose = args["verbose"];

        std::cerr << "Clustering consensus sequences..." << std::endl;
        auto clusters = cluster_reads(reads, 6, 0.5, 25, 0.4, 0.4, 0.05, 0, false, 0.15, is_rna, verbose, n_threads, false);
        auto correction = correct_reads(clusters, reads, 0.3, 0.3, 30.0, 200, 0, n_threads, verbose);

        int cid = 0;
        for (auto &r: correction.consensi) {
            int total_reads = 0;
            auto creads = clusters[cid].seqs;

            for (auto &s: creads) {
                auto info = split(reads[s.seq_id].header, '=');
                int rcount = std::stoi(info[1]);
                total_reads += rcount;
            }

            r.header += " total_reads=" + std::to_string(total_reads);
            cid++;
        }

        write_fastq_file(correction.consensi, args["output"].as<std::string>(".") + "/transcriptome.fq");
    } else {
        std::cerr << "Unknown mode. More info" << std::endl;
    }
    
}