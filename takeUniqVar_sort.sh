#!/bin/sh

# ============================================================================
# Shell script to get unique variants from their variant annotation
# ============================================================================

############### VARIABLES TO CHANGE ###############
INPUT_FILE=$1
OUTPUT_FILE=$2

DELIM=$'\t'
VAR_COL=11
###################################################

#Keep the header and write to ouput file
head -1 $INPUT_FILE > $OUTPUT_FILE

#Sort the rest of the body and append to output file
tail -n +2 $INPUT_FILE | sort -u -t"$DELIM" -k"$VAR_COL","$VAR_COL"  >> $OUTPUT_FILE
