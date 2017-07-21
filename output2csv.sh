#!/bin/sh

# =============================================================================
# Script that takes in a tab-delimited .output file, converts all pre-existing
# ',' to ';', then converts all tabs to commas (',')
# =============================================================================

INPUT_FILE=$1
OUTPUT_FILE=$2


less $INPUT_FILE | tr ',' ';' | tr $'\t' ',' > $OUTPUT_FILE
