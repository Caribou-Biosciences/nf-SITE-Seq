// Copyright 2025 Caribou Biosciences. All rights reserved.

#include <unistd.h>

#include <algorithm>
#include <deque>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <queue>
#include <string>
#include <vector>

#include "common.h"

/*
This program reads in distributions of read counts and computes a sliding window
difference statistic over all positions in the genome. The input distributions
are in TSV format, are potentially bgzip-compressed, and contain read counts
from one or more sample replicates at each position. Sample replicates are
paired so the output distribution has the same number of measurements at each
position as the input distributions. The program outputs the resulting
distribution in TSV format to stdout and only includes records for which at
least one measurement is non-zero. The program also collapses adjacent positions
with the same stat values into a single record.

Normally the window difference statistic is calculated from two separate input
distributions. However, the selfComparison CLI parameter enables a single input
distribution to be used in which pairs of sample replicates within the same
distribution are used as the inputs to the statistic.

The program takes advantage of the fact that the inputs are sorted by coordinate
to read each input record once.
*/

struct ComboPositionRecord {
    // This struct tracks a combination of test and control values at a single
    // position
    std::string chromosome;
    int64_t position;
    std::vector<double> testCounts;
    std::vector<double> controlCounts;
};

// Apply rotation resampling on vectors, i.e. make the 2nd element become the
// first, etc
template <typename T>
void applyRotationResampling(std::vector<T> *vector) {
    if ((*vector).size() > 1) {
        std::rotate((*vector).begin(), (*vector).begin() + 1, (*vector).end());
    }
}

// Add or remove counts from running sums
void updateSums(const ComboPositionRecord &comboPos,
                std::vector<double> *testSums, std::vector<double> *controlSums,
                bool isAdding) {
    for (size_t idx = 0; idx < comboPos.testCounts.size(); idx++) {
        if (isAdding) {
            (*testSums)[idx] += comboPos.testCounts[idx];
            (*controlSums)[idx] += comboPos.controlCounts[idx];
        } else {
            (*testSums)[idx] -= comboPos.testCounts[idx];
            (*controlSums)[idx] -= comboPos.controlCounts[idx];
        }
    }
}

// Calculate window difference statistic from test and control sum vectors
void calculateWindowDiffs(std::vector<double> *windowDiffs,
                          const std::vector<double> &testSums,
                          const std::vector<double> &controlSums) {
    for (size_t idx = 0; idx < (*windowDiffs).size(); idx++) {
        (*windowDiffs)[idx] = std::max(testSums[idx] - controlSums[idx], 0.0);
    }
}

void emitRecord(std::ostringstream &buffer, const std::string &chrom,
                int64_t start, int64_t end, const std::vector<double> &vals,
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
    for (const double &val : vals) {
        double roundedVal = roundToDecimalPlaces(val, precision);
        buffer << "\t" << roundedVal;
    }
    buffer << "\n";

    // Flush buffer if it reaches the limit
    if (buffer.tellp() >= BUFFER_LIMIT) {
        std::cout << buffer.str();
        buffer.str("");
        buffer.clear();
    }
}

void ingestRecordIntoHeap(
    std::priority_queue<PositionRecord, std::vector<PositionRecord>,
                        std::greater<>> *minHeap,
    const RangeRecord &record, const bool isControl, const bool rotate) {
    // Ingest a record into the heap. Input records cover a genomic range, we
    // expand the record into each individual position within the range and push
    // it into the heap. If requested, apply rotation resampling to make the
    // second element the first element
    PositionRecord pos = {record.chromosome, record.start, record.counts,
                          isControl, false};
    if (rotate) {
        applyRotationResampling(&pos.counts);
    }
    for (int64_t i = record.start; i < record.end; i++) {
        pos.position = i;
        pos.isLastPos = (i == record.end - 1);
        (*minHeap).push(pos);
    }
}

void printUsage(const char *programName) {
    std::cerr << "Usage: " << programName
              << " -t <test.tsv[.gz]> -c <control.tsv[.gz]> -w <W> -f <faidx>\n"
              << "  -t <file>   Test distribution file (required)\n"
              << "  -c <file>   Control distribution file (required unless -s "
                 "is used)\n"
              << "  -w <int>    Window size (required)\n"
              << "  -f <file>   FASTA index file (required)\n"
              << "  -s          Use self-comparison mode\n"
              << "  -p <int>    Decimal precision for output (default: 6)\n";
}

int main(int argc, char *argv[]) {
    // Parse options
    std::string testFile;
    std::string controlFile;
    std::string faidxFile;
    int W = -1;
    bool selfComparison = false;
    int precision = DEFAULT_PRECISION;
    int opt;
    while ((opt = getopt(argc, argv, "t:c:w:f:sp::")) != -1) {
        switch (opt) {
            case 't':
                testFile = optarg;
                break;
            case 'c':
                controlFile = optarg;
                break;
            case 'w':
                W = std::stoi(optarg);
                break;
            case 'f':
                faidxFile = optarg;
                break;
            case 's':
                selfComparison = true;
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
    if (testFile.empty() || faidxFile.empty() || W <= 0 ||
        (controlFile.empty() && !selfComparison)) {
        std::cerr << "Error: Parameters -t, -f, and -w must be provided, and "
                     "-c is required unless -s is provided\n";
        printUsage(argv[0]);
        return 1;
    }

    // Open input distribution files, get sample names
    std::unique_ptr<DistTSVFile> testReader =
        std::make_unique<DistTSVFile>(testFile);
    int numSamples = testReader->sampleNames.size();
    std::unique_ptr<DistTSVFile> controlReader;
    std::vector<std::string> testSamples, controlSamples;
    if (!selfComparison) {
        // Not a self comparison, load the control file
        controlReader = std::make_unique<DistTSVFile>(controlFile);
        if (controlReader->sampleNames.size() >
            testReader->sampleNames.size()) {
            std::cerr << "Warning: Control file has more columns than test "
                         "file, extra control columns will be ignored\n";
        } else if (controlReader->sampleNames.size() !=
                   testReader->sampleNames.size()) {
            throw std::invalid_argument(
                "Test and control files have different numbers of samples");
        }
        testSamples = testReader->sampleNames;
        controlSamples = controlReader->sampleNames;
    } else {
        // If self comparison, ensure enough test samples and assign control
        // headers as rotated test headers
        testSamples = testReader->sampleNames;
        if (testSamples.size() < 2) {
            std::cerr << "Error: Input file must have at least 2 samples when "
                         "performing a self-comparison\n";
            return 1;
        }
        controlSamples.assign(testSamples.begin(), testSamples.end());
        applyRotationResampling(&controlSamples);
    }

    // Create output headers
    std::vector<std::string> outputSampleHeaders;
    outputSampleHeaders.reserve(numSamples);
    for (size_t i = 0; i < numSamples; ++i) {
        outputSampleHeaders.push_back(testSamples[i] + "~" + controlSamples[i]);
    }

    // Get contig sizes from faidx
    std::map<std::string, int64_t> contigSizes = parseFaidx(faidxFile);

    // Setup output buffer
    std::ostringstream buffer;
    buffer << std::fixed << std::setprecision(precision);

    // Write header to output
    buffer << "#chromosome\tstart\tend";
    for (const auto &header : outputSampleHeaders) {
        buffer << "\t" << header;
    }
    buffer << "\n";

    // Initialize data structures for processing
    // Use a heap to store positions and easily get the one with the minimum
    // coordinate
    std::priority_queue<PositionRecord, std::vector<PositionRecord>,
                        std::greater<>>
        minHeap;
    std::deque<ComboPositionRecord> window;
    std::vector<double> runningTestSums(numSamples,
                                        0.0);  // current test window sums
    std::vector<double> runningControlSums(numSamples,
                                           0.0);  // current control window sums
    std::vector<double> runningWindowDiffs(numSamples,
                                           0.0);  // current window diffs
    std::vector<double> lastWindowDiffs(
        numSamples, 0.0);  // prev window diffs, is output next

    RangeRecord testRecord, controlRecord;
    PositionRecord testPos, controlPos;
    ComboPositionRecord curComboPos;

    std::string curChrom;
    std::string prevChrom = "nowaythiscouldbearealcontignameright?";
    int64_t currentStart = -1;  // start position of the next record to output

    // Read initial records and populate heap
    bool testAvailable = testReader->readRecord(&testRecord);
    if (testAvailable) {
        ingestRecordIntoHeap(&minHeap, testRecord, false, false);
    }

    bool controlAvailable;
    if (selfComparison) {
        controlAvailable = testAvailable;
        if (controlAvailable) {
            ingestRecordIntoHeap(&minHeap, testRecord, true,
                                 true);  // ingest with rotation resampling
        }
    } else {
        controlAvailable = controlReader->readRecord(&controlRecord);
        if (controlAvailable) {
            ingestRecordIntoHeap(&minHeap, controlRecord, true, false);
        }
    }

    // Main processing loop
    while (!minHeap.empty()) {
        // Get the minimum position from the heap
        PositionRecord minPos = minHeap.top();
        minHeap.pop();

        curChrom = minPos.chromosome;
        if (curChrom != prevChrom) {
            auto entry = contigSizes.find(curChrom);
            if (entry == contigSizes.end()) {
                throw std::invalid_argument("Contig not found in faidx: " +
                                            curChrom);
            }
        }

        // Determine if this is test or control record
        if (minPos.isControl) {
            controlPos = minPos;
        } else {
            testPos = minPos;
        }

        // Check if there's a matching position in the heap
        if (!minHeap.empty() && minHeap.top().chromosome == minPos.chromosome &&
            minHeap.top().position == minPos.position) {
            PositionRecord matchPos = minHeap.top();
            minHeap.pop();
            if (matchPos.isControl == minPos.isControl) {
                throw std::runtime_error(
                    "Duplicate record found for position " +
                    std::to_string(minPos.position) + " in " +
                    std::string(minPos.isControl ? "control" : "test"));
            } else if (matchPos.isControl) {
                controlPos = matchPos;
            } else {
                testPos = matchPos;
            }

            // Read next record if needed
            if (matchPos.isLastPos) {
                if (!selfComparison) {
                    if (matchPos.isControl && controlAvailable) {
                        controlAvailable =
                            controlReader->readRecord(&controlRecord);
                        if (controlAvailable) {
                            ingestRecordIntoHeap(&minHeap, controlRecord, true,
                                                 false);
                        }
                    } else if (!matchPos.isControl && testAvailable) {
                        testAvailable = testReader->readRecord(&testRecord);
                        if (testAvailable) {
                            ingestRecordIntoHeap(&minHeap, testRecord, false,
                                                 false);
                        }
                    }
                }
            }
        } else {
            // No matching position for other sample, assume 0 counts
            if (minPos.isControl) {
                testPos = {minPos.chromosome, minPos.position,
                           std::vector<double>(numSamples, 0.0), false};
            } else {
                controlPos = {minPos.chromosome, minPos.position,
                              std::vector<double>(numSamples, 0.0), true};
            }
        }

        // Combine test and control data
        curComboPos = {testPos.chromosome, testPos.position, testPos.counts,
                       controlPos.counts};

        // Remove positions that will be outside the window when curComboPos is
        // added
        while (!window.empty() &&
               (window.front().chromosome != curComboPos.chromosome ||
                window.front().position < curComboPos.position - 2 * W)) {
            const ComboPositionRecord toRemove = window.front();
            window.pop_front();

            // Update sums and calculate new window diffs
            updateSums(toRemove, &runningTestSums, &runningControlSums, false);
            calculateWindowDiffs(&runningWindowDiffs, runningTestSums,
                                 runningControlSums);

            // Emit record if window diffs have changed
            if (runningWindowDiffs != lastWindowDiffs &&
                (toRemove.chromosome != curComboPos.chromosome ||
                 toRemove.position < curComboPos.position - 2 * W - 1)) {
                int64_t endPos = toRemove.position + W + 1;
                emitRecord(buffer, toRemove.chromosome, currentStart, endPos,
                           lastWindowDiffs, precision,
                           contigSizes[toRemove.chromosome]);
                // Reset window diffs and record start position
                lastWindowDiffs = runningWindowDiffs;
                currentStart = endPos;
            }
        }

        // Add current position to window
        window.push_back(curComboPos);

        if (curChrom != prevChrom) {
            // Reset for new chromosome
            currentStart = curComboPos.position - W;
            runningTestSums = curComboPos.testCounts;
            runningControlSums = curComboPos.controlCounts;
            calculateWindowDiffs(&runningWindowDiffs, runningTestSums,
                                 runningControlSums);
            lastWindowDiffs = runningWindowDiffs;
        } else {
            // Update running sums and calculate window diffs
            updateSums(curComboPos, &runningTestSums, &runningControlSums,
                       true);
            calculateWindowDiffs(&runningWindowDiffs, runningTestSums,
                                 runningControlSums);
        }

        // Emit record if window diffs have changed
        if (runningWindowDiffs != lastWindowDiffs) {
            int64_t endPos = curComboPos.position - W;
            emitRecord(buffer, curComboPos.chromosome, currentStart, endPos,
                       lastWindowDiffs, precision,
                       contigSizes[curComboPos.chromosome]);
            lastWindowDiffs = runningWindowDiffs;
            currentStart = endPos;
        }

        // Read next record if needed
        if (minPos.isLastPos) {
            if (selfComparison && testAvailable) {
                testAvailable = testReader->readRecord(&testRecord);
                controlAvailable = testAvailable;
                if (testAvailable) {
                    ingestRecordIntoHeap(&minHeap, testRecord, false, false);
                    ingestRecordIntoHeap(&minHeap, testRecord, true, true);
                }
            } else if (minPos.isControl && controlAvailable) {
                controlAvailable = controlReader->readRecord(&controlRecord);
                if (controlAvailable) {
                    ingestRecordIntoHeap(&minHeap, controlRecord, true, false);
                }
            } else if (!minPos.isControl && testAvailable) {
                testAvailable = testReader->readRecord(&testRecord);
                if (testAvailable) {
                    ingestRecordIntoHeap(&minHeap, testRecord, false, false);
                }
            }
        }

        prevChrom = curChrom;
    }

    // Process any remaining positions in the window
    while (!window.empty()) {
        const ComboPositionRecord toRemove = window.front();
        window.pop_front();

        updateSums(toRemove, &runningTestSums, &runningControlSums, false);
        calculateWindowDiffs(&runningWindowDiffs, runningTestSums,
                             runningControlSums);

        if (runningWindowDiffs != lastWindowDiffs) {
            int64_t endPos = toRemove.position + W + 1;
            emitRecord(buffer, toRemove.chromosome, currentStart, endPos,
                       lastWindowDiffs, precision,
                       contigSizes[toRemove.chromosome]);
            lastWindowDiffs = runningWindowDiffs;
            currentStart = endPos;
        }
    }

    // Output any remaining data in buffer
    std::cout << buffer.str();
    return 0;
}
