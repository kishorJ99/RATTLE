#include "fasta.hpp"
#include <fstream>
#include <algorithm>
#include <iostream>
#include <zlib.h>

std::string unzip_file(std::string filename, int index){
    gzFile infile = gzopen(filename.c_str(), "rb");
    std::cerr << "Start decompressing file" << std::endl;
    filename = filename.substr(0, index);
    FILE *outfile = fopen(filename.c_str(), "wb");
    gzrewind(infile);

    if (infile == NULL) {
        std::cerr << "Error: Failed to decompress the file" << std::endl;
        return " ";
    }

    unsigned char unzipBuffer[1048576];
    std::vector<unsigned char> unzippedData;
    while(!gzeof(infile))
    {
        int len = gzread(infile, unzipBuffer, 1048576);
        fwrite(unzipBuffer, 1, len, outfile);
    }
    fclose(outfile);
    gzclose(infile);

    std::cerr << "Decompressing file Complete" << std::endl;
    return filename;
}

read_set_t read_fasta_file(std::string file, bool raw, int lower_len, int upper_len) {
    read_set_t result;

    std::ifstream infile(file);
    std::string line;
    std::string header;
    std::string seq;
    std::string qt = "";
    bool isLinux = true;

    // Setting the FASTA format quality value to the maximum value '~' for correction and polish step
    std::getline(infile, line);
    if(char (line[line.size() - 1]) == '\r'){
        isLinux = false;
        header = line.substr(0, line.size() - 1);
    } else {
        header = line;
    }

    if(isLinux){
        while (std::getline(infile, line)) {
            if (line.size() == 0) continue;

            if (line[0] == '>') {
                if (!header.empty()) {
                    std::transform(seq.begin(), seq.end(),seq.begin(), ::toupper);
                    if(raw){
                        for(int i = 0; i < seq.size(); i++){
                            qt += '~';
                        }
                        read_t r{header, seq, "", qt};
                        result.push_back(r);
                        qt = "";
                    } else if( seq.length() >= lower_len && seq.length() <= upper_len){
                        for(int i = 0; i < seq.size(); i++){
                            qt += '~';
                        }
                        read_t r{header, seq, "", qt};
                        result.push_back(r);
                        qt = "";
                    } 
                }

                seq = "";
                header = line;
            } else {
                seq += line;
            }
        }
    } else {
        while (std::getline(infile, line)) {
            if (line.size() == 0) continue;

            if (line[0] == '>') {
                if (!header.empty()) {
                    std::transform(seq.begin(), seq.end(),seq.begin(), ::toupper);
                    if(raw){
                        for(int i = 0; i < seq.size(); i++){
                            qt += '~';
                        }
                        read_t r{header, seq, "", qt};
                        result.push_back(r);
                        qt = "";
                    } else if( seq.length() >= lower_len && seq.length() <= upper_len){
                        for(int i = 0; i < seq.size(); i++){
                            qt += '~';
                        }
                        read_t r{header, seq, "", qt};
                        result.push_back(r);
                        qt = "";
                    } 
                }

                seq = "";
                header = line.substr(0, line.size() - 1);
            } else {
                seq += line.substr(0, line.size() - 1);
            }
         }     
    }

    std::transform(seq.begin(), seq.end(),seq.begin(), ::toupper);
    if(raw){
        for(int i = 0; i < seq.size(); i++){
            qt += '~';
        }
        read_t r{header, seq, "", qt};
        result.push_back(r);
        qt = "";
    } else if( seq.length() >= lower_len && seq.length() <= upper_len){
        for(int i = 0; i < seq.size(); i++){
            qt += '~';
        }
        read_t r{header, seq, "", qt};
        result.push_back(r);
        qt = "";
    } 
    
    return result;
}

read_set_t read_fastq_file(std::string file, bool raw, int lower_len, int upper_len) {
    read_set_t result;

    std::ifstream infile(file);
    std::string line;
    std::string header;
    std::string seq;
    std::string ann;
    std::string qt;
    int lineID = 0;
    bool isLinux = true;

    std::getline(infile, line);
    if(line[line.size() - 1] == '\r'){
        isLinux = false;
        header = line.substr(0, line.size() - 1);
        lineID++;
    } else {
        header = line;
        lineID++;
    }

    if(isLinux){
        while (std::getline(infile, line)) {
            if (lineID == 0) {
                header = line;
                ++lineID;
            } else if (lineID == 1) {
                seq = line;
                ++lineID;
            } else if (lineID == 2) {
                ann = line;
                ++lineID;
            } else if (lineID == 3) {
                qt = line;
                lineID = 0;

                if(raw){
                    read_t r{header, seq, ann, qt};
                    result.push_back(r);
                } else if( seq.length() >= lower_len && seq.length() <= upper_len){
                    read_t r{header, seq, ann, qt};
                    result.push_back(r);
                }
            }
        }
    } else {
        while (std::getline(infile, line)) {
            if (lineID == 0) {
                header = line.substr(0, line.size() - 1);
                ++lineID;
            } else if (lineID == 1) {
                seq = line.substr(0, line.size() - 1);
                ++lineID;
            } else if (lineID == 2) {
                ann = line.substr(0, line.size() - 1);
                ++lineID;
            } else if (lineID == 3) {
                qt = line.substr(0, line.size() - 1);
                lineID = 0;

                if(raw){
                    read_t r{header, seq, ann, qt};
                    result.push_back(r);
                } else if( seq.length() >= lower_len && seq.length() <= upper_len){
                    read_t r{header, seq, ann, qt};
                    result.push_back(r);
                }
            }
        }
    }
    return result;
}

void write_fastq_file(const read_set_t &reads, std::string file) {
    std::ofstream f;
    f.open(file);

    for (auto &r: reads) {
        f << r.header << std::endl << r.seq << std::endl << r.ann << std::endl << r.quality << std::endl;
    }

    f.close();
}

bool _comp_read_set_desc(read_t a, read_t b) {
    return a.seq.size() > b.seq.size();
}

void sort_read_set(read_set_t &rs) {
    std::stable_sort(rs.begin(), rs.end(), _comp_read_set_desc);
}