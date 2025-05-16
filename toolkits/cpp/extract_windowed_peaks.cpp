// Copyright 2025 Caribou Biosciences. All rights reserved.

#include <unistd.h>

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "common.h"

/*
This program reads in a distribution of values across a reference genome and
identifies peaks of the distribution. The input distribution is in TSV format,
is potentially bgzip-compressed, and contains values from one or more sample
replicates at each position. The peak identification process is very simple - a
sliding window is passed over the distribution and records are partitioned into
groups where each record in a group is within W base pairs of at least one other
record in the group. The midpoint of the record with the maximum sum of values
is chosen as the peak location and is written to stdout in TSV format. The
program also creates a file summarizing how many times each total and meean sum
peak value is observed; this is output as a TSV file.

The program takes advantage of the fact that the input is sorted by coordinate
to read each record once.
*/

void emitRecord(std::ostringstream &buffer, const std::string chrom,
                const int64_t start, const int64_t end,
                const std::vector<double> &vals) {
    buffer << chrom << "\t" << start << "\t" << end;
    for (size_t i = 0; i < vals.size(); ++i) {
        buffer << "\t" << vals[i];
    }
    buffer << "\n";

    // Flush buffer if it reaches the limit
    if (buffer.tellp() >= BUFFER_LIMIT) {
        std::cout << buffer.str();
        buffer.str("");
        buffer.clear();
    }
}

void findPeaks(const std::string &inputFile, const int windowSize,
               const std::string &outputDistributionFile, const int precision) {
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

    // Write header to output
    buffer << "#chromosome\tstart\tend";
    for (const auto &sampleName : inputReader->sampleNames) {
        buffer << "\t" << sampleName;
    }
    buffer << "\n";

    // Track how many times each peak value is observed
    std::unordered_map<double, int64_t> peakDistribution;

    // Read first record
    RangeRecord curRecord;
    bool recordAvailable = inputReader->readRecord(&curRecord);

    std::string curChrom;    // Track current chromosome
    int64_t lastPos = -1;    // Track last position observed
    double curMaxSum = -1;   // Track the maximum value in the current window
    RangeRecord peakRecord;  // Track the record with the maximum value in the
                             // current window
    while (recordAvailable) {
        if (curRecord.chromosome != curChrom ||
            curRecord.start > lastPos + windowSize) {
            // Current record is outside of window from last position observed
            if (lastPos != -1) {
                // Output current max. Set peak location to center of max
                // record.
                int64_t midPos =
                    peakRecord.start + (peakRecord.end - peakRecord.start) / 2;
                emitRecord(buffer, peakRecord.chromosome, midPos, midPos + 1,
                           peakRecord.counts);
                peakDistribution[roundToDecimalPlaces(curMaxSum, precision)]++;
            }
            // Reset window
            curMaxSum = curRecord.totalCount();
            peakRecord = curRecord;
        } else {
            // Current record is within window from last position observed.
            // Check if it is the new maximum
            if (curRecord.totalCount() >= curMaxSum) {
                peakRecord = curRecord;
                curMaxSum = curRecord.totalCount();
            }
        }

        if (curRecord.totalCount() > 0.0) {
            curChrom = curRecord.chromosome;
            lastPos = curRecord.end;
        }

        // Read next record
        recordAvailable = inputReader->readRecord(&curRecord);
    }

    // Handle last window
    if (lastPos != -1) {
        int64_t midPos =
            peakRecord.start + (peakRecord.end - peakRecord.start) / 2;
        emitRecord(buffer, peakRecord.chromosome, midPos, midPos + 1,
                   peakRecord.counts);
        peakDistribution[roundToDecimalPlaces(curMaxSum, precision)]++;
    }

    std::cout << buffer.str();

    // Write the distribution of peak sums/means. Order by descending sum.
    std::vector<std::pair<double, int64_t>> sortedDistribution;
    std::copy(peakDistribution.begin(), peakDistribution.end(),
              std::back_inserter(sortedDistribution));
    auto comparator = [](const std::pair<double, int64_t> &a,
                         const std::pair<double, int64_t> &b) {
        return a.first > b.first;
    };
    std::sort(sortedDistribution.begin(), sortedDistribution.end(), comparator);
    distributionOutput << "sum\tmean\tcount" << '\n';
    for (const auto &entry : sortedDistribution) {
        const double mean = entry.first / numSamples;
        distributionOutput << entry.first << '\t' << std::fixed
                           << std::setprecision(6) << mean << '\t'
                           << entry.second << '\n';
    }
    distributionOutput.close();
}

void printUsage(const char *programName) {
    std::cerr << "Usage: " << programName
              << " -i <input.tsv[.gz]> -d output.dist.tsv -w <W>\n"
              << "  -i <file>   Input file (required)\n"
              << "  -d <file>   Output distribution file (required)\n"
              << "  -w <int>    Window size (required)\n"
              << "  -p <int>    Decimal precision for output (default: 6)\n";
}

int main(int argc, char *argv[]) {
    std::string inputFile, outputDistributionFile;
    int windowSize = -1;
    int precision = DEFAULT_PRECISION;

    int opt;
    while ((opt = getopt(argc, argv, "i:w:d:p::")) != -1) {
        switch (opt) {
            case 'i':
                inputFile = optarg;
                break;
            case 'w':
                windowSize = std::stoi(optarg);
                break;
            case 'd':
                outputDistributionFile = optarg;
                break;
            case 'p':
                precision = std::stoi(optarg);
                break;
            default:
                printUsage(argv[0]);
                return 1;
        }
    }

    // Validate required parameters
    if (inputFile.empty() || outputDistributionFile.empty() || windowSize < 0) {
        std::cerr << "Parameters -i, -w, and -d are required\n";
        printUsage(argv[0]);
        return 1;
    }

    findPeaks(inputFile, windowSize, outputDistributionFile, precision);

    return 0;
}
