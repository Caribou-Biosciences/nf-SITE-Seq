// Copyright 2025 Caribou Biosciences. All rights reserved.

#include <unistd.h>

#include <algorithm>
#include <deque>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "common.h"

/*
This program reads in a distribution of read counts and computes a sliding
window sum over all positions in the genome. The input distribution is in TSV
format, is potentially bgzip-compressed, and contains read counts from one or
more sample replicates at each position. The window sum is calculated on each
replicate so the output distribution has the same number of measurements at each
position as the input distribution. The program outputs the resulting
distribution in TSV format to stdout and only includes records for which at
least one measurement is non-zero. The program also collapses adjacent positions
with the same stat values into a single record. The program also creates a file
summarizing how many times each total and meean sum is observed; this is output
as a TSV file.

The program takes advantage of the fact that the input is sorted by coordinate
to read each record once.
*/

void emitRecord(std::ostringstream &buffer, const std::string &chrom,
                int64_t start, int64_t end, const std::vector<double> &vals,
                std::unordered_map<double, int64_t> *sumDistribution,
                const int precision, const int64_t contigSize) {
    // Only emit records with a valid contig position and at least one non-zero
    // value. If coordinates are partially out of bounds, clip to bounds.
    if (start < 0) {
        if (end > 0) {
            start = 0;
        } else {
            return;
        }
    }
    if (end > contigSize) {
        if (start < contigSize) {
            end = contigSize;
        } else {
            return;
        }
    }
    double total = std::accumulate(vals.begin(), vals.end(), 0.0);
    if (areAlmostEqual(total, 0.0, precision)) {
        return;
    }

    buffer << chrom << "\t" << start << "\t" << std::min(end, contigSize);
    for (size_t i = 0; i < vals.size(); ++i) {
        double rounded = roundToDecimalPlaces(vals[i], precision);
        buffer << "\t" << rounded;
    }
    buffer << "\n";

    // Add total to distribution counts
    double roundedTotal = roundToDecimalPlaces(total, precision);
    (*sumDistribution)[roundedTotal] += end - start;

    // Flush buffer if it reaches the limit
    if (buffer.tellp() >= BUFFER_LIMIT) {
        std::cout << buffer.str();
        buffer.str("");
        buffer.clear();
    }
}

void printUsage(const char *programName) {
    std::cerr << "Usage: " << programName
              << " -i <input.tsv[.gz]> -d output.dist.tsv -w <W> -f <faidx>\n"
              << "  -i <file>   Input file (required)\n"
              << "  -d <file>   Output distribution file (required)\n"
              << "  -w <int>    Window size (required)\n"
              << "  -f <file>   FASTA index file (required)\n"
              << "  -p <int>    Decimal precision for output (default: 6)\n";
}

int main(int argc, char *argv[]) {
    // Parse options
    std::string inputFile;
    std::string outputDistributionFile;
    std::string faidxFile;
    int W = -1;
    int precision = DEFAULT_PRECISION;
    int opt;
    while ((opt = getopt(argc, argv, "i:f:w:d:p::")) != -1) {
        switch (opt) {
            case 'i':
                inputFile = optarg;
                break;
            case 'w':
                W = std::stoi(optarg);
                break;
            case 'd':
                outputDistributionFile = optarg;
                break;
            case 'p':
                precision = std::stoi(optarg);
                break;
            case 'f':
                faidxFile = optarg;
                break;
            default:
                printUsage(argv[0]);
                return 1;
        }
    }

    // Validate required parameters
    if (inputFile.empty() || outputDistributionFile.empty() ||
        faidxFile.empty() || W <= 0) {
        std::cerr << "Error: Parameters -i, -f, -d and -w must be provided\n";
        printUsage(argv[0]);
        return 1;
    }

    // Open input file
    std::unique_ptr<DistTSVFile> inputReader =
        std::make_unique<DistTSVFile>(inputFile);
    int numSamples = inputReader->sampleNames.size();

    // Open output file for distribution counts
    std::ofstream distributionOutput(outputDistributionFile);
    distributionOutput << std::fixed << std::setprecision(precision);

    // Setup output buffer
    std::ostringstream buffer;
    buffer << std::fixed << std::setprecision(precision);

    // Get contig sizes from faidx
    std::map<std::string, int64_t> contigSizes = parseFaidx(faidxFile);

    // Write header to output
    buffer << "#chromosome\tstart\tend";
    for (const auto &sampleName : inputReader->sampleNames) {
        buffer << "\t" << sampleName;
    }
    buffer << "\n";

    // Initialize data structures for processing
    std::deque<PositionRecord> window;  // positions in the current window
    std::vector<double> runningSums(numSamples, 0.0);  // current window sums
    std::vector<double> lastSums(numSamples,
                                 0.0);  // prev window sums, is output next
    std::unordered_map<double, int64_t>
        sumDistribution;  // map storing how many times each sum is observed

    // Record tracking structs
    RangeRecord curRangeRec;
    PositionRecord curPosRec;

    // Track current and previous chromosomes
    std::string curChrom;
    std::string prevChrom = "nowaythiscouldbearealcontignameright?";
    int64_t currentStart = -1;  // start position of the next record to output

    // Read first record
    bool recordAvailable = inputReader->readRecord(&curRangeRec);

    // Main processing loop
    while (recordAvailable) {
        curChrom = curRangeRec.chromosome;
        if (curChrom != prevChrom) {
            auto entry = contigSizes.find(curChrom);
            if (entry == contigSizes.end()) {
                throw std::invalid_argument("Contig not found in faidx: " +
                                            curChrom);
            }
        }

        // Iterate through each position in the current range record
        curPosRec = {curRangeRec.chromosome, curRangeRec.start,
                     curRangeRec.counts};
        for (int64_t position = curRangeRec.start; position < curRangeRec.end;
             position++) {
            curPosRec.position = position;

            // Remove positions that will be outside the window when curPosRec
            // is added
            while (!window.empty() &&
                   (window.front().chromosome != curPosRec.chromosome ||
                    window.front().position < curPosRec.position - 2 * W)) {
                const PositionRecord toRemove = window.front();
                window.pop_front();

                // Remove position from running sums
                for (size_t i = 0; i < numSamples; i++) {
                    runningSums[i] -= toRemove.counts[i];
                }

                // Emit record if window sums have changed
                if (runningSums != lastSums &&
                    (toRemove.chromosome != curPosRec.chromosome ||
                     toRemove.position < curPosRec.position - 2 * W - 1)) {
                    int64_t endPos = toRemove.position + W + 1;
                    emitRecord(buffer, toRemove.chromosome, currentStart,
                               endPos, lastSums, &sumDistribution, precision,
                               contigSizes[toRemove.chromosome]);
                    // Reset window sums and record start position
                    lastSums = runningSums;
                    currentStart = endPos;
                }
            }

            // Add current position to window
            window.push_back(curPosRec);
            if (curChrom != prevChrom) {
                // Reset for new chromosome
                currentStart = curPosRec.position - W;
                runningSums = curPosRec.counts;
                lastSums = curPosRec.counts;
            } else {
                // Same chromsome, update running window sums
                for (size_t i = 0; i < numSamples; i++) {
                    runningSums[i] += curPosRec.counts[i];
                }
            }

            // Emit record if window sums have changed
            if (runningSums != lastSums) {
                int64_t endPos = curPosRec.position - W;
                emitRecord(buffer, curPosRec.chromosome, currentStart, endPos,
                           lastSums, &sumDistribution, precision,
                           contigSizes[curPosRec.chromosome]);
                // Reset window sums and record start position
                lastSums = runningSums;
                currentStart = endPos;
            }

            prevChrom = curChrom;
        }

        // Read a new record if available
        if (recordAvailable)
            recordAvailable = inputReader->readRecord(&curRangeRec);
    }

    // Process any remaining positions in the window
    while (!window.empty()) {
        const PositionRecord toRemove = window.front();
        for (size_t i = 0; i < numSamples; i++) {
            runningSums[i] -= toRemove.counts[i];
        }
        window.pop_front();

        if (runningSums != lastSums) {
            int64_t endPos = toRemove.position + W + 1;
            emitRecord(buffer, toRemove.chromosome, currentStart, endPos,
                       lastSums, &sumDistribution, precision,
                       contigSizes[toRemove.chromosome]);
            lastSums = runningSums;
            currentStart = endPos;
        }
    }

    std::cout << buffer.str();

    // Write the distribution of peak sums/means. Order by descending sum.
    std::vector<std::pair<double, int64_t>> sortedDistribution;
    std::copy(sumDistribution.begin(), sumDistribution.end(),
              std::back_inserter(sortedDistribution));
    auto comparator = [](const std::pair<double, int64_t> &a,
                         const std::pair<double, int64_t> &b) {
        return a.first > b.first;
    };
    std::sort(sortedDistribution.begin(), sortedDistribution.end(), comparator);
    distributionOutput << "sum\tmean\tcount" << '\n';
    for (const auto &entry : sortedDistribution) {
        const double mean = static_cast<double>(entry.first) / numSamples;
        distributionOutput << entry.first << '\t' << mean << '\t'
                           << entry.second << '\n';
    }
    distributionOutput.close();

    return 0;
}
