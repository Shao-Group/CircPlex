#!/bin/bash

# Check if input file is provided
if [[ $# -ne 1 ]]; then
    echo "Usage: $0 <fasta_file>"
    exit 1
fi

INPUT_FILE="$1"
TEMP_FILE="${INPUT_FILE}.tmp"

# Check if input file exists
if [[ ! -f "$INPUT_FILE" ]]; then
    echo "Error: File '$INPUT_FILE' not found!"
    exit 1
fi

# Update sequence length in headers, handling "norepeat" and "nocycle" cases
awk '
/^>/ {
    if (seq != "") {
        if (seq == "norepeat" || seq == "nocycle" || seq == "toolong") {
            print header " unitlength=0"
        } else {
            print header " unitlength=" length(seq)
        }
        print seq
        seq = ""
    }
    header = $0
    next
}
{
    seq = seq $0
}
END {
    if (seq != "") {
        if (seq == "norepeat" || seq == "nocycle" || seq == "toolong") {
            print header " unitlength=0"
        } else {
            print header " unitlength=" length(seq)
        }
        print seq
    }
}' "$INPUT_FILE" > "$TEMP_FILE"

# Replace the original file with the updated one
mv "$TEMP_FILE" "$INPUT_FILE"

echo "Sequence lengths updated in $INPUT_FILE"
