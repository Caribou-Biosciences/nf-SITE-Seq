// Copyright 2025 Caribou Biosciences. All rights reserved.

#include <htslib/bgzf.h>
#include <htslib/hts.h>
#include <unistd.h>

#include <algorithm>
#include <fstream>
#include <iostream>
#include <memory>
#include <nlohmann/json.hpp>
#include <string>
#include <vector>

#include "common.h"

/*
This program reads in two paired .fastq.gz files and filters reads based on
the quality scores of their UMI bases. The UMI is assumed to be at the start
of read 2 with a length defined by an input parameter. A read is filtered out
if it has more than N bases below a quality score threshold Q, where N and Q
are defined by the input parameters.
The program also creates a report in JSON format containing the number of
input/output/failing reads.
*/

// Function to determine if a read passes the quality filter
bool passes_quality_filter(const std::string &qualities, const int umi_length,
                           const int phred_threshold,
                           const int max_bases_below) {
    if (qualities.empty() || umi_length <= 0) {
        return false;
    }

    int below_threshold_count = 0;
    const int check_length =
        std::min(static_cast<size_t>(umi_length), qualities.size());

    for (int i = 0; i < check_length; ++i) {
        if (qualities[i] - '!' < phred_threshold) {
            if (++below_threshold_count > max_bases_below) {
                return false;
            }
        }
    }
    return true;
}

// Helper function for writing buffered data to BGZF files
bool write_to_bgzf(BGZF *file, std::ostringstream &buffer) {
    if (!file) {
        std::cerr << "Error: BGZF file pointer is null\n";
        return false;
    }

    if (buffer.tellp() >= static_cast<std::streampos>(BUFFER_LIMIT)) {
        std::string bufferData = buffer.str();
        if (bgzf_write(file, bufferData.c_str(), bufferData.size()) < 0) {
            std::cerr << "Error writing to BGZF file\n";
            return false;
        }
        buffer.str("");
        buffer.clear();
    }
    return true;
}

// Parse command line arguments
struct ProgramOptions {
    std::string read1_path;
    std::string read2_path;
    std::string output_read1_path;
    std::string output_read2_path;
    std::string output_report_path;
    int umi_length = 0;
    int phred_threshold = 0;
    int max_bases_below_threshold = 0;

    bool validate() const {
        return !read1_path.empty() && !read2_path.empty() &&
               !output_read1_path.empty() && !output_read2_path.empty() &&
               !output_report_path.empty() && umi_length > 0 &&
               phred_threshold > 0 && max_bases_below_threshold >= 0;
    }

    void print_usage(const char *program_name) const {
        std::cerr << "Usage: " << program_name
                  << " -1 <read1.fastq[.gz]> -2 <read2.fastq[.gz]>"
                  << " -o <output_read1.fastq.gz> -O <output_read2.fastq.gz>"
                  << " -r <output_report.json> -U <UMI_length>"
                  << " -P <PHRED_threshold> -M <max_bases_below_threshold>\n";
    }
};

ProgramOptions parse_args(int argc, char *argv[]) {
    ProgramOptions options;
    int opt;

    while ((opt = getopt(argc, argv, "1:2:o:O:r:U:P:M:")) != -1) {
        switch (opt) {
            case '1':
                options.read1_path = optarg;
                break;
            case '2':
                options.read2_path = optarg;
                break;
            case 'o':
                options.output_read1_path = optarg;
                break;
            case 'O':
                options.output_read2_path = optarg;
                break;
            case 'r':
                options.output_report_path = optarg;
                break;
            case 'U':
                options.umi_length = std::stoi(optarg);
                break;
            case 'P':
                options.phred_threshold = std::stoi(optarg);
                break;
            case 'M':
                options.max_bases_below_threshold = std::stoi(optarg);
                break;
            default:
                options.print_usage(argv[0]);
                exit(1);
        }
    }

    if (!options.validate()) {
        std::cerr << "Error: Missing or invalid arguments.\n";
        options.print_usage(argv[0]);
        exit(1);
    }

    return options;
}

// Custom RAII class for htsFile pointers
class HtsFileHandle {
 private:
    htsFile *file;

 public:
    explicit HtsFileHandle(const std::string &path, const char *mode)
        : file(hts_open(path.c_str(), mode)) {}
    ~HtsFileHandle() {
        if (file) hts_close(file);
    }
    htsFile *get() { return file; }
    bool is_valid() const { return file != nullptr; }
    HtsFileHandle(const HtsFileHandle &) = delete;
    HtsFileHandle &operator=(const HtsFileHandle &) = delete;
};

// Custom RAII class for BGZF pointers
class BgzfHandle {
 private:
    BGZF *file;

 public:
    explicit BgzfHandle(const std::string &path, const char *mode)
        : file(bgzf_open(path.c_str(), mode)) {}
    ~BgzfHandle() {
        if (file) bgzf_close(file);
    }
    BGZF *get() { return file; }
    bool is_valid() const { return file != nullptr; }
    BgzfHandle(const BgzfHandle &) = delete;
    BgzfHandle &operator=(const BgzfHandle &) = delete;
};

// Custom RAII class for kstring_t
class KStringHandle {
 private:
    kstring_t str;

 public:
    KStringHandle() : str({0, 0, nullptr}) {}
    ~KStringHandle() { free(str.s); }
    kstring_t *get() { return &str; }
    KStringHandle(const KStringHandle &) = delete;
    KStringHandle &operator=(const KStringHandle &) = delete;
};

int main(int argc, char *argv[]) {
    auto options = parse_args(argc, argv);

    // Open input files with RAII
    HtsFileHandle read1_file(options.read1_path, "r");
    HtsFileHandle read2_file(options.read2_path, "r");

    if (!read1_file.is_valid() || !read2_file.is_valid()) {
        std::cerr << "Error: Unable to open input files\n";
        return 1;
    }

    // Open output files with RAII
    BgzfHandle output_read1_file(options.output_read1_path, "w");
    BgzfHandle output_read2_file(options.output_read2_path, "w");

    if (!output_read1_file.is_valid() || !output_read2_file.is_valid()) {
        std::cerr << "Error: Unable to open output files\n";
        return 1;
    }

    // Counters
    int total_reads = 0;
    int passing_reads = 0;
    int failing_reads = 0;

    // Use RAII for kstring_t
    KStringHandle read1_buffer;
    KStringHandle read2_buffer;

    // Buffers for bulk writing
    std::ostringstream write_buffer_read1, write_buffer_read2;

    // Process reads
    while (hts_getline(read1_file.get(), '\n', read1_buffer.get()) >= 0 &&
           hts_getline(read2_file.get(), '\n', read2_buffer.get()) >= 0) {
        // Parse read 1
        std::string read1_header = read1_buffer.get()->s;
        hts_getline(read1_file.get(), '\n', read1_buffer.get());
        std::string read1_sequence = read1_buffer.get()->s;
        hts_getline(read1_file.get(), '\n', read1_buffer.get());
        std::string read1_plus = read1_buffer.get()->s;
        hts_getline(read1_file.get(), '\n', read1_buffer.get());
        std::string read1_quality = read1_buffer.get()->s;

        // Parse read 2
        std::string read2_header = read2_buffer.get()->s;
        hts_getline(read2_file.get(), '\n', read2_buffer.get());
        std::string read2_sequence = read2_buffer.get()->s;
        hts_getline(read2_file.get(), '\n', read2_buffer.get());
        std::string read2_plus = read2_buffer.get()->s;
        hts_getline(read2_file.get(), '\n', read2_buffer.get());
        std::string read2_quality = read2_buffer.get()->s;

        total_reads++;

        // Filter reads based on quality
        if (passes_quality_filter(read2_quality, options.umi_length,
                                  options.phred_threshold,
                                  options.max_bases_below_threshold)) {
            passing_reads++;
            write_buffer_read1 << read1_header << '\n'
                               << read1_sequence << '\n'
                               << read1_plus << '\n'
                               << read1_quality << '\n';
            write_buffer_read2 << read2_header << '\n'
                               << read2_sequence << '\n'
                               << read2_plus << '\n'
                               << read2_quality << '\n';

            if (!write_to_bgzf(output_read1_file.get(), write_buffer_read1) ||
                !write_to_bgzf(output_read2_file.get(), write_buffer_read2)) {
                return 1;
            }
        } else {
            failing_reads++;
        }
    }

    // Final write of any remaining buffered data
    const std::string &read1_remaining = write_buffer_read1.str();
    if (!read1_remaining.empty()) {
        if (bgzf_write(output_read1_file.get(), read1_remaining.c_str(),
                       read1_remaining.size()) < 0) {
            std::cerr
                << "Error writing remaining data to BGZF file for read1\n";
            return 1;
        }
    }

    const std::string &read2_remaining = write_buffer_read2.str();
    if (!read2_remaining.empty()) {
        if (bgzf_write(output_read2_file.get(), read2_remaining.c_str(),
                       read2_remaining.size()) < 0) {
            std::cerr
                << "Error writing remaining data to BGZF file for read2\n";
            return 1;
        }
    }

    // Write JSON report
    using json = nlohmann::json;
    json report = {
        {"total_reads", total_reads},
        {"passing_reads", passing_reads},
        {"failing_reads", failing_reads},
        {"passing_percentage",
         total_reads > 0 ? (100.0 * passing_reads / total_reads) : 0.0}};

    std::ofstream json_file(options.output_report_path);
    if (!json_file) {
        std::cerr << "Error: Unable to open report file for writing\n";
        return 1;
    }
    json_file << report.dump(4);  // Pretty-print with 4 spaces

    return 0;
}
