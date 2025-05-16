# Use an official lightweight base image
FROM ubuntu:20.04

# Set environment variables for non-interactive installs
ENV DEBIAN_FRONTEND=noninteractive

# Install dependencies
RUN apt-get update && \
    apt-get install -y \
    build-essential \
    wget \
    curl \
    gdb \
    valgrind \
    libbz2-dev \
    liblzma-dev \
    zlib1g-dev \
    libcurl4-gnutls-dev \
    software-properties-common \
    && rm -rf /var/lib/apt/lists/*

RUN add-apt-repository ppa:deadsnakes/ppa && \
    apt-get update && \
    apt-get install -y \
    python3.10 \
    python3.10-venv \
    python3.10-dev

RUN rm /usr/bin/python3 && \
    ln -s python3.10 /usr/bin/python3

RUN curl -sS https://bootstrap.pypa.io/get-pip.py | python3.10

RUN pip3 install --upgrade pip

# Install HTSlib
RUN wget https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2 && \
    tar -xvjf htslib-1.16.tar.bz2 && \
    cd htslib-1.16 && \
    ./configure && \
    make && \
    make install && \
    cd .. && rm -rf htslib-1.16*

# Set LD_LIBRARY_PATH to include HTSlib
ENV LD_LIBRARY_PATH="/usr/local/lib"

# Download json library
RUN mkdir -p /usr/local/include/nlohmann
RUN curl -L -o /usr/local/include/nlohmann/json.hpp https://github.com/nlohmann/json/raw/refs/tags/v3.11.3/single_include/nlohmann/json.hpp

COPY python /usr/src/siteseq_python
RUN pip3 install /usr/src/siteseq_python

# Copy the C++ source code into the container
COPY cpp /usr/src/siteseq_cpp

# Compile the C++ program
RUN g++ -std=c++17 -o /usr/local/bin/merge_read_starts /usr/src/siteseq_cpp/merge_read_starts.cpp /usr/src/siteseq_cpp/common.cpp -lhts
RUN g++ -std=c++17 -o /usr/local/bin/extract_windowed_peaks /usr/src/siteseq_cpp/extract_windowed_peaks.cpp /usr/src/siteseq_cpp/common.cpp -lhts
RUN g++ -std=c++17 -o /usr/local/bin/calculate_window_diffs /usr/src/siteseq_cpp/calculate_window_diffs.cpp /usr/src/siteseq_cpp/common.cpp -lhts
RUN g++ -std=c++17 -o /usr/local/bin/quality_filter_umis /usr/src/siteseq_cpp/quality_filter_umis.cpp /usr/src/siteseq_cpp/common.cpp -lhts
RUN g++ -std=c++17 -o /usr/local/bin/calculate_window_sums /usr/src/siteseq_cpp/calculate_window_sums.cpp /usr/src/siteseq_cpp/common.cpp -lhts
RUN g++ -std=c++17 -o /usr/local/bin/count_read_starts /usr/src/siteseq_cpp/count_read_starts.cpp /usr/src/siteseq_cpp/common.cpp -lhts
