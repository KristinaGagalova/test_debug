#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <reads_1> <reads_2>"
    exit 1
fi

R1="$1"
R2="$2"


if [ ! -f "$R1" ]; then
    echo "Error: File '$R1' not found!"
    exit 1
fi

if [ ! -f "$R2" ]; then
    echo "Error: File '$R2' not found!"
    exit 1
fi



if [[ "$R1" =~ \.gz$ ]]; then
    echo "Processing gzipped file: $R1"
    zcat "$R1" | sed 's:\(.*\).* \(.*\)\( length=.*\):\1-\2/1:g' | gzip > "$R1.tmp" && mv "$R1.tmp" "$R1"
else
    echo "Processing uncompressed file: $R1"
    sed -i 's:\(.*\).* \(.*\)\( length=.*\):\1-\2/1:g' "$R1"
fi
echo "Processed: $R1"


if [[ "$R2" =~ \.gz$ ]]; then
    echo "Processing gzipped file: $R2"
    zcat "$R2" | sed 's:\(.*\).* \(.*\)\( length=.*\):\1-\2/2:g' | gzip > "$R2.tmp" && mv "$R2.tmp" "$R2"
else
    echo "Processing uncompressed file: $R2"
    sed -i 's:\(.*\).* \(.*\)\( length=.*\):\1-\2/2:g' "$R2"
fi
echo "Processed: $R2"

echo "All files processed successfully!"