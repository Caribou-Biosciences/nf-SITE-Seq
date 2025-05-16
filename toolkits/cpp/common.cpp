// Copyright 2025 Caribou Biosciences. All rights reserved.

#include "common.h"

#include <map>
#include <string>

/*
This file stores common variables, structs, classes, and functions
*/

// Global variables
const size_t BUFFER_LIMIT = 2000000;
const size_t DEFAULT_PRECISION = 10;

// Inline function definitions for RangeRecord struct

double RangeRecord::totalCount() const {
    return std::accumulate(counts.begin(), counts.end(), 0.0);
}

// Inline function definitions for PositionRecord struct

bool PositionRecord::operator>(const PositionRecord &other) const {
    // Allows PositionRecord struct to be sorted by coordinate
    return std::tie(chromosome, position) >
           std::tie(other.chromosome, other.position);
}

double PositionRecord::totalCount() const {
    return std::accumulate(counts.begin(), counts.end(), 0.0);
}

// Utility functions

bool isBgzipped(const std::string &filename) {
    return (filename.size() > 3 &&
            filename.substr(filename.size() - 3) == ".gz") ||
           (filename.size() > 4 &&
            filename.substr(filename.size() - 4) == ".bgz");
}

bool getLineFromBGZF(BGZF *fp, std::string *line) {
    // Get a line from a bgzipped file
    kstring_t str = {0, 0, NULL};
    int ret = bgzf_getline(fp, '\n', &str);
    if (ret < 0) {
        ks_free(&str);
        return false;
    }
    (*line).assign(str.s);
    ks_free(&str);
    return true;
}

std::string getBaseNameFromPath(const std::string &path) {
    // Get the basename of a path, i.e. the filename with all extensions removed
    size_t lastSlash = path.find_last_of("/\\");
    std::string basename =
        (lastSlash == std::string::npos) ? path : path.substr(lastSlash + 1);
    size_t firstDot = basename.find_first_of('.');
    if (firstDot != std::string::npos) {
        basename.resize(firstDot);
    }
    return basename;
}

std::map<std::string, int64_t> parseFaidx(const std::string &filename) {
    // Parse a faidx file to get a mapping of contig names to sizes
    std::map<std::string, int64_t> contigSizes;
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file: " + filename);
    }
    std::string line;
    while (std::getline(file, line)) {
        std::istringstream ss(line);
        std::string contigName;
        int64_t contigSize;
        if (ss >> contigName >> contigSize) {
            contigSizes[contigName] = contigSize;
        }
    }
    return contigSizes;
}

bool areAlmostEqual(double a, double b, int precision) {
    double epsilon = std::pow(10.0, -precision);
    return std::fabs(a - b) < epsilon;
}

double roundToDecimalPlaces(const double value, const int precision) {
    // Round the input value to a specified number of decimal places and
    // ensure negative zero is avoided
    double factor = std::pow(10.0, precision);
    double rounded = std::round(value * factor) / factor;

    // Avoid negative zero
    if (rounded == 0) {
        rounded = 0.0;  // Ensure it is positive zero
    }
    return rounded;
}

// Inline function definitions for DistTSVFile class

DistTSVFile::DistTSVFile(const std::string &file) : filename(file) {
    baseName = getBaseNameFromPath(filename);
    isBgzipCompressed = isBgzipped(file);
    if (isBgzipCompressed) {
        bgzfFile = bgzf_open(file.c_str(), "r");
        if (!bgzfFile) {
            throw std::invalid_argument("Could not open bgzip file: " + file);
        }
    } else {
        fileStream.open(file);
        if (!fileStream.is_open()) {
            throw std::invalid_argument("Could not open file: " + file);
        }
    }
    getSampleNames(&sampleNames);
}

void DistTSVFile::getSampleNames(std::vector<std::string> *names) {
    std::string line;
    bool success = false;
    if (isBgzipCompressed) {
        success = getLineFromBGZF(bgzfFile, &line);
    } else {
        success = static_cast<bool>(std::getline(fileStream, line));
    }
    if (success) {
        std::stringstream ss(line);
        std::string header;
        std::getline(ss, header, '\t');  // chrom
        std::getline(ss, header, '\t');  // start
        std::getline(ss, header, '\t');  // end
        while (std::getline(ss, header, '\t')) {
            (*names).push_back(header);  // sample name
        }
    } else {
        throw std::invalid_argument("Could not parse file");
    }
}

bool DistTSVFile::readRecord(RangeRecord *record) {
    std::string line;
    bool success = false;
    if (isBgzipCompressed) {
        success = getLineFromBGZF(bgzfFile, &line);
    } else {
        success = static_cast<bool>(std::getline(fileStream, line));
    }
    if (success) {
        if (line.empty() ||
            line.find_first_not_of(" \t\n\v\f\r") == std::string::npos) {
            return false;
        }
        (*record).counts.clear();
        std::istringstream iss(line);
        if (!(iss >> (*record).chromosome >> (*record).start >>
              (*record).end)) {
            return false;
        }
        double count;
        while (iss >> count) {
            (*record).counts.push_back(count);
        }
        if ((*record).counts.size() < sampleNames.size()) {
            throw std::invalid_argument(
                "Incorrect number of sample columns in line");
        }
        return true;
    }
    return false;
}

DistTSVFile::~DistTSVFile() {
    if (bgzfFile) {
        bgzf_close(bgzfFile);
    }
}
