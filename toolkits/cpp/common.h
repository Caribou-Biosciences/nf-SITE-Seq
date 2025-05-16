// Copyright 2025 Caribou Biosciences. All rights reserved.

#ifndef TOOLKITS_CPP_COMMON_H_
#define TOOLKITS_CPP_COMMON_H_

#include <htslib/bgzf.h>
#include <htslib/kstring.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

/*
This file stores common variables, structs, classes, and functions
*/

// Global variables
extern const size_t BUFFER_LIMIT;
extern const size_t DEFAULT_PRECISION;

// Record structs

struct RangeRecord {
    // This struct tracks input records - a set of values defined over a
    // coordinate range
    std::string chromosome;
    int64_t start;
    int64_t end;
    std::vector<double> counts;

    double totalCount() const;
};

struct PositionRecord {
    // This struct tracks a single input position
    std::string chromosome;
    int64_t position;
    std::vector<double> counts;
    bool isControl;
    bool isLastPos;
    int index;

    bool operator>(const PositionRecord &other) const;

    double totalCount() const;
};

// Utility functions

bool isBgzipped(const std::string &filename);

bool getLineFromBGZF(BGZF *fp, std::string *line);

std::string getBaseNameFromPath(const std::string &path);

std::map<std::string, int64_t> parseFaidx(const std::string &filename);

bool areAlmostEqual(double a, double b, int precision);

double roundToDecimalPlaces(const double value, const int precision);

// Class for parsing distribution TSVs
class DistTSVFile {
 private:
    std::ifstream fileStream;
    BGZF *bgzfFile = nullptr;
    bool isBgzipCompressed;
    std::string filename;

    void getSampleNames(std::vector<std::string> *names);

 public:
    std::string baseName;
    std::vector<std::string> sampleNames;

    explicit DistTSVFile(const std::string &file);

    bool readRecord(RangeRecord *record);

    ~DistTSVFile();
};

#endif  // TOOLKITS_CPP_COMMON_H_
