#!/bin/bash
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
PROGRAM_NAME="site-seq-tools"
VERSION="1.0.0"
docker build $SCRIPT_DIR -t "$PROGRAM_NAME:$VERSION"
