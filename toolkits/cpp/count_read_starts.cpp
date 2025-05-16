// Copyright 2025 Caribou Biosciences. All rights reserved.

#include <getopt.h>
#include <htslib/sam.h>

#include <iostream>
#include <map>
#include <nlohmann/json.hpp>
#include <string>

#include "common.h"

/*
This program parses a BAM file and counts the number of reads that start at each
position in the genome. The program writes the results to stdout in TSV format
with "chromosome", "start", "end", and "count" columns. Several filters are
applied to the reads before they are included in the output distribution: (1)
minimum reference alignment length, (2) minimum mapping quality, (3) maximum
clipping at the start of the alignment. If requested, output counts are
normalized by the total number of aligned reads per million. Finally, a report
in JSON format is created with the number of input/output reads and the number
of reads passing each filter.
*/

void printUsage(const char *programName) {
    std::cerr
        << "Usage: " << programName
        << " -b <input.bam> -s <output.stats.json>\n"
        << "  -b <file>   Input BAM file (required)\n"
        << "  -s <file>   Output summary stats JSON file path (required)\n"
        << "  -l <int>    Minimum alignment length (default: 0)\n"
        << "  -c <int>    Maximum clipping length (default: 0)\n"
        << "  -c <int>    Minimum mapping quality (default: 0)\n"
        << "  -n          Apply normalization (default: false)\n"
        << "  -p <int>    Decimal precision for output (default: 6)\n";
}

int main(int argc, char *argv[]) {
    char *bamPath = nullptr;
    char *statsPath = nullptr;
    int minAlignmentLen = 0;
    int maxClipping = 0;
    int minMapQual = 0;
    bool applyNorm = false;
    int precision = DEFAULT_PRECISION;

    struct option longOptions[] = {
        {"bam-path", required_argument, nullptr, 'b'},
        {"stats-path", required_argument, nullptr, 's'},
        {"min-aln-length", required_argument, nullptr, 'l'},
        {"max-clipping", required_argument, nullptr, 'c'},
        {"min-map-qual", required_argument, nullptr, 'q'},
        {"apply-norm", optional_argument, nullptr, 'n'},
        {"precision", optional_argument, nullptr, 'p'},
        {nullptr, 0, nullptr, 0}};

    int opt;
    int optionIndex = 0;
    while ((opt = getopt_long(argc, argv, "b:s:l:c:q:np::", longOptions,
                              &optionIndex)) != -1) {
        try {
            switch (opt) {
                case 'b':
                    bamPath = optarg;
                    break;
                case 's':
                    statsPath = optarg;
                    break;
                case 'l':
                    minAlignmentLen = std::stoi(optarg);
                    break;
                case 'c':
                    maxClipping = std::stoi(optarg);
                    break;
                case 'q':
                    minMapQual = std::stoi(optarg);
                    break;
                case 'n':
                    applyNorm = true;
                    break;
                case 'p':
                    precision = std::stoi(optarg);
                    break;
                case '?':
                    std::cerr << "Invalid option or missing argument\n";
                    return 1;
                default:
                    printUsage(argv[0]);
                    break;
            }
        } catch (const std::exception &e) {
            std::cerr << "Error parsing command-line arguments: " << e.what()
                      << "\n";
            printUsage(argv[0]);
            return 1;
        }
    }

    // Validate required input paths
    if (!bamPath || !statsPath) {
        std::cerr << "Error: --bam-path and --stats-path are required.\n";
        printUsage(argv[0]);
        return 1;
    }

    // Open BAM file
    samFile *in = sam_open(bamPath, "r");
    if (!in) {
        std::cerr << "Error: Could not open BAM file.\n";
        return 1;
    }

    // Initialize report output file
    std::ofstream statsFile(statsPath);
    if (!statsFile) {
        std::cerr << "Error: Unable to open output stats file." << std::endl;
        sam_close(in);
        return 1;
    }

    // Load BAM header and initialize HTSlib structures
    bam_hdr_t *header = sam_hdr_read(in);
    bam1_t *read = bam_init1();

    // A map to store read counts at each position per chromosome
    std::map<std::string, std::map<int64_t, int64_t>> positionCounts;

    // Variables to track read counts
    int64_t totalReads = 0, numReadsPassFilters = 0;
    int64_t numBelowMapQ = 0, numAboveClipping = 0, numBelowRefLength = 0;

    // Process each read in BAM file
    while (sam_read1(in, header, read) >= 0) {
        int32_t tid = read->core.tid;  // Chromosome ID
        int64_t pos = read->core.pos;  // Start position of the alignment
        bool is_reverse =
            read->core.flag & BAM_FREVERSE;  // Strand of alignment

        if (tid < 0 || pos < 0) continue;

        // Make sure this is a primary alignment
        uint16_t flag = read->core.flag;
        if (flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue;

        totalReads++;

        // Parse CIGAR string to (1) determine how long the alignment is on the
        // reference and (2) determine how much clipping there is at the start
        // of the read
        uint32_t *cigar = bam_get_cigar(read);
        int n_cigar = read->core.n_cigar;
        int startClip = 0, refLength = 0;
        for (int i = 0; i < n_cigar; ++i) {
            int op = bam_cigar_op(cigar[i]);
            int len = bam_cigar_oplen(cigar[i]);

            // Check clipping. If forward strand, first operation is at start of
            // read; if reverse strand, last operation is at start of read
            if (i == 0 && (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) &&
                !is_reverse) {
                startClip = len;
            } else if (i == n_cigar - 1 &&
                       (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) &&
                       is_reverse) {
                startClip = len;
            }

            // Increment reference length counter when CIGAR operation consumes
            // reference bases
            if (op == BAM_CMATCH || op == BAM_CDEL || op == BAM_CREF_SKIP ||
                op == BAM_CEQUAL || op == BAM_CDIFF) {
                refLength += len;
            }
        }

        // Check filters
        bool passesFilters = true;
        if (read->core.qual < minMapQual) {
            numBelowMapQ++;
            passesFilters = false;
        }
        if (startClip > maxClipping) {
            numAboveClipping++;
            passesFilters = false;
        }
        if (refLength < minAlignmentLen) {
            numBelowRefLength++;
            passesFilters = false;
        }

        if (!passesFilters) continue;
        numReadsPassFilters++;

        // Get coordinate of start position and increment count at that
        // position. If reverse strand, start position is the end of the
        // alignment
        std::string chrom = header->target_name[tid];
        int64_t start_pos = is_reverse ? pos + refLength - 1 : pos;
        positionCounts[chrom][start_pos]++;
    }

    // Close BAM file and free memory
    bam_destroy1(read);
    bam_hdr_destroy(header);
    sam_close(in);

    // If requested, normalize by per-million reads
    double normFactor;
    if (applyNorm) {
        normFactor = static_cast<double>(totalReads) / 1000000;
    } else {
        normFactor = 1.0;
    }

    // Buffer for accumulating output
    std::ostringstream buffer;
    buffer << std::fixed << std::setprecision(precision);
    // Write TSV header
    buffer << "#chromosome\tstart\tend\tcount\n";

    // Write each position with its read count to stdout
    int64_t totalSites = 0;
    for (const auto &chromEntry : positionCounts) {
        const std::string &chrom = chromEntry.first;
        for (const auto &posEntry : chromEntry.second) {
            totalSites++;
            int64_t pos = posEntry.first;
            double normCount =
                static_cast<double>(posEntry.second) / normFactor;
            // Output half-open intervals
            buffer << chrom << "\t" << pos << "\t" << pos + 1 << "\t"
                   << normCount << "\n";

            // Flush buffer if it reaches the limit
            if (buffer.tellp() >= BUFFER_LIMIT) {
                std::cout << buffer.str();
                buffer.str("");
                buffer.clear();
            }
        }
    }
    std::cout << buffer.str();

    // Create JSON report
    using json = nlohmann::json;
    json statSummary = {{"total_reads", totalReads},
                        {"passing_reads", numReadsPassFilters},
                        {"total_sites", totalSites},
                        {"filters",
                         {{"low_mapping_quality", numBelowMapQ},
                          {"above_clipping", numAboveClipping},
                          {"below_ref_length", numBelowRefLength}}}};
    statsFile << statSummary.dump(4);  // Pretty-print with 4 spaces
    statsFile.close();

    return 0;
}
