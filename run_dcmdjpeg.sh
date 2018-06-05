#!/bin/bash
# define input and output dirs
# INPUT_DIR="/some/dir/where/the/dicoms/are/stored"
INPUT_DIR=$1
OUTPUT_DIR="$INPUT_DIR/decompressed_dicoms"
# create output dir
mkdir -p "$OUTPUT_DIR"
# call dcmdjpeg with each file in input dir
#for i in $(find "$INPUT_DIR" -name '*.dcm'); do
for i in $(find "$INPUT_DIR" -type f); do	
# if the dicom files do not have the extension .dcm use instead: for i in $(find "$INPUT_DIR" -type f); do
 FILENAME="${i#$INPUT_DIR/}"
 #FILEPARENT="${FILENAME%/*}"
 #if [[ "$FILEPARENT" != "$FILENAME" -a -N "$FILEPARENT" ]]; then
 # # create new parent folders if needed
 # mkdir -p "$FILEPARENT"
 #fi
 dcmdjpeg "$i" "$OUTPUT_DIR/$FILENAME"
done