// Copyright 2025 Caribou Biosciences. All rights reserved.

#include <getopt.h>

#include <algorithm>
#include <functional>
#include <iomanip>
#include <iostream>
#include <memory>
#include <queue>
#include <string>
#include <unordered_set>
#include <vector>

#include "common.h"

/*
This program reads in multiple files containing distributions of values across a
reference genome combines them into a single distribution that is output to
stdout. The input distributions are in TSV format, are bgzip-compressed, and
contain values from a single sample replicate at each position. The program
outputs a single record for each input position containing the values from each
of the input samples at that position, filling with zeros when a sample has no
value at a position.

The program takes advantage of the fact that the input is sorted by coordinate
to read each record once.
*/

void printUsage(const char *programName) {
    std::cerr << "Usage: " << programName
              << " -b <input1.bam> <input2.bam> ... <inputN.bam>\n"
              << "  -b <file1> <file2> ... <fileN>  Input files (required)\n"
              << "  -p <int>                        Decimal precision for "
                 "output (default: 6)\n";
}

int main(int argc, char *argv[]) {
    std::vector<std::string> tsvPaths;
    int precision = DEFAULT_PRECISION;
    struct option longOptions[] = {
        {"bam-paths", required_argument, nullptr, 'b'},
        {"precision", optional_argument, nullptr, 'p'},
        {nullptr, 0, nullptr, 0}};

    int opt;
    int optionIndex = 0;
    while ((opt = getopt_long(argc, argv, "b:p::", longOptions,
                              &optionIndex)) != -1) {
        switch (opt) {
            case 'b':
                for (int i = optind - 1; i < argc; ++i) {
                    if (argv[i][0] == '-') {
                        break;
                    }
                    tsvPaths.push_back(argv[i]);
                    ++optind;
                }
                optind = std::min(optind, argc - 1);
                break;

            case 'p':
                precision = std::stoi(optarg);
                break;

            case '?':
                printUsage(argv[0]);
                return 1;

            default:
                printUsage(argv[0]);
                break;
        }
    }

    if (tsvPaths.size() < 1) {
        std::cerr << "At least one BAM file must be provided" << std::endl;
        printUsage(argv[0]);
        return 1;
    }

    std::unordered_set<std::string> nameSet;  // Used to check duplicate names
    std::vector<std::unique_ptr<DistTSVFile>> readers(
        tsvPaths.size());  // Used to store file readers

    // Record variables
    RangeRecord rangeRec;
    PositionRecord positionRec, otherPositionRec, nextPositionRec;

    // Use a heap to store positions and easily get the one with the minimum
    // coordinate
    std::priority_queue<PositionRecord, std::vector<PositionRecord>,
                        std::greater<>>
        minHeap;

    // Go through each input file, create file reader, get file base name, read
    // the first data record, and ingest single positions from the first record
    // into the heap
    for (int i = 0; i < tsvPaths.size(); ++i) {
        std::string tsvPath = tsvPaths[i];
        readers[i] = std::make_unique<DistTSVFile>(tsvPath);
        if (nameSet.find(readers[i]->baseName) != nameSet.end()) {
            std::cerr << "Error: Name collision detected for '"
                      << readers[i]->baseName << "'\n";
            return 1;
        } else if (readers[i]->sampleNames.size() != 1) {
            std::cerr << "Input TSV does not have a single sample column: '"
                      << tsvPath << "'\n";
            return 1;
        }
        nameSet.insert(readers[i]->baseName);

        if (readers[i]->readRecord(&rangeRec)) {
            positionRec = {rangeRec.chromosome,
                           rangeRec.start,
                           rangeRec.counts,
                           false,
                           false,
                           i};
            for (int64_t pos = rangeRec.start; pos < rangeRec.end; pos++) {
                positionRec.position = pos;
                if (pos == rangeRec.end - 1) {
                    positionRec.isLastPos = true;
                }
                minHeap.push(positionRec);
            }
        }
    }

    // Setup output buffer
    std::ostringstream buffer;
    buffer << std::fixed << std::setprecision(precision);

    // Write header to output
    buffer << "#chromosome\tstart\tend";
    for (const auto &reader : readers) {
        buffer << "\t" << reader->baseName;
    }
    buffer << "\n";

    while (!minHeap.empty()) {
        // Get the position record at the top of the heap - this will be the
        // next output position
        positionRec = minHeap.top();
        minHeap.pop();

        // Initialize a vector storing the value at this position from each
        // input file
        std::vector<double> counts(tsvPaths.size(), 0.0);
        counts[positionRec.index] = positionRec.counts[0];

        // Get any position records at the same position from the heap
        while (!minHeap.empty() &&
               minHeap.top().chromosome == positionRec.chromosome &&
               minHeap.top().position == positionRec.position) {
            otherPositionRec = minHeap.top();
            minHeap.pop();
            counts[otherPositionRec.index] =
                otherPositionRec.counts[0];  // Add the count for this input

            // If this was the last position in a range record, read a new
            // record from the same input and ingest into heap
            if (otherPositionRec.isLastPos &&
                readers[otherPositionRec.index]->readRecord(&rangeRec)) {
                nextPositionRec = {rangeRec.chromosome,
                                   rangeRec.start,
                                   rangeRec.counts,
                                   false,
                                   false,
                                   otherPositionRec.index};
                for (int64_t pos = rangeRec.start; pos < rangeRec.end; pos++) {
                    nextPositionRec.position = pos;
                    if (pos == rangeRec.end - 1) {
                        nextPositionRec.isLastPos = true;
                    }
                    minHeap.push(nextPositionRec);
                }
            }
        }

        // Output current position
        buffer << positionRec.chromosome << "\t" << positionRec.position << "\t"
               << positionRec.position + 1;
        for (const auto &count : counts) {
            buffer << "\t" << count;
        }
        buffer << "\n";

        // If this was the last position in a range record, read a new record
        // from the same input and ingest into heap
        if (positionRec.isLastPos &&
            readers[positionRec.index]->readRecord(&rangeRec)) {
            nextPositionRec = {rangeRec.chromosome,
                               rangeRec.start,
                               rangeRec.counts,
                               false,
                               false,
                               positionRec.index};
            for (int64_t pos = rangeRec.start; pos < rangeRec.end; pos++) {
                nextPositionRec.position = pos;
                if (pos == rangeRec.end - 1) {
                    nextPositionRec.isLastPos = true;
                }
                minHeap.push(nextPositionRec);
            }
        }

        // Flush buffer if it reaches the limit
        if (buffer.tellp() >= BUFFER_LIMIT) {
            std::cout << buffer.str();
            buffer.str("");
            buffer.clear();
        }
    }

    std::cout << buffer.str();

    return 0;
}
