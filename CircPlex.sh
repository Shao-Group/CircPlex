#!/bin/bash
# ============================================================
# CircPlex: A circular RNA assembler from complex rolling circular long reads
# ============================================================
# Usage:
# ./CircPlex.sh <input.fasta/fastq> <outprefix> <kmc-path> <equirep-path> <minimap2-path>
# Example:
# ./CircPlex.sh ./example/input.fasta ./example/output ./KMC/bin ./equirep-1.0.0/src ./minimap2-2.24_x64-linux 
# ============================================================

# ----------------------------
# Arguments
# ----------------------------
if [ "$#" -ne 5 ]; then
    echo "Usage: $0 <input.fasta/fastq> <outprefix> <kmc-path> <equirep-path> <minimap2-path>"
    exit 1
fi

INPUT="$1"
OUTPREFIX="$2"
KMC_DIR="$3"
EQUIREP_DIR="$4"
MINIMAP2_DIR="$5"

# ----------------------------
# Config
# ----------------------------
SCRIPT_DIR="./src"
K=10

# ----------------------------
# Functions
# ----------------------------
function info() {
    echo -e "\033[1;34m[INFO]\033[0m $1"
}

function error() {
    echo -e "\033[1;31m[ERROR]\033[0m $1" >&2
    exit 1
}

function check_file() {
    if [ ! -f "$1" ]; then
        error "File not found: $1"
    fi
}

# Check mandatory files
check_file "$INPUT"
check_file "$KMC_DIR/kmc"
check_file "$KMC_DIR/kmc_tools"
check_file "$EQUIREP_DIR/EquiRep"
check_file "$MINIMAP2_DIR/minimap2"

# ----------------------------
# Prepare working directory
# ----------------------------
INPUT_BASENAME=$(basename "$INPUT")
INPUT_NAME="${INPUT_BASENAME%.*}"
WORKING_DIR="./${INPUT_NAME}_tmp"
mkdir -p "$WORKING_DIR"

info "Input file: $INPUT"
info "Working directory: $WORKING_DIR"
info "Output prefix: $OUTPREFIX"

# ----------------------------
# Step 1: Compile C++ programs
# ----------------------------
info "[1/9] Compiling C++ programs..."
g++ -O3 -ffast-math -march=native "$SCRIPT_DIR/find_A.cpp" -o "$SCRIPT_DIR/find_A"
g++ -O3 -ffast-math -march=native "$SCRIPT_DIR/make_two_copy.cpp" -o "$SCRIPT_DIR/make_two_copy"
g++ -O3 -ffast-math -march=native "$SCRIPT_DIR/calculate_cigar.cpp" -o "$SCRIPT_DIR/calculate_cigar"
g++ -O3 -ffast-math -march=native "$SCRIPT_DIR/combine.cpp" -o "$SCRIPT_DIR/combine"
g++ -O3 -ffast-math -march=native "$SCRIPT_DIR/check_complexity.cpp" -o "$SCRIPT_DIR/check_complexity"

# ----------------------------
# Step 2: Convert FASTA to FASTQ if needed
# ----------------------------
FASTQ_FILE="$WORKING_DIR/${INPUT_NAME}.fastq"
if [[ "$INPUT" == *.fasta || "$INPUT" == *.fa ]]; then
    info "[2/9] Converting FASTA to FASTQ..."
    awk 'BEGIN{OFS=""} /^>/ {print "@"substr($0,2); getline seq; print seq; print "+"; print gensub(/./, "I", "g", seq)}' "$INPUT" > "$FASTQ_FILE"
else
    info "[2/9] Input is FASTQ, copying to working directory..."
    cp "$INPUT" "$FASTQ_FILE"
fi

# ----------------------------
# Step 3: Run per-read KMC and detect complex reads
# ----------------------------
info "[3/9] Running per-read KMC counting and complexity detection..."

KMC_COUNTS_DIR="${WORKING_DIR}/kmer_counts"
mkdir -p "$KMC_COUNTS_DIR"
COMPLEX_FASTA="$WORKING_DIR/${INPUT_NAME}_complex_reads.fasta"
> "$COMPLEX_FASTA"

readnum=0
while read -r header && read -r seq && read -r plus && read -r qual; do
    ((readnum++))
    readname="${header#@}"
    readname="${readname%% *}"
    header_extra="${header#@${readname} }"

    TMP_R_FASTQ="${KMC_COUNTS_DIR}/tmp_${readnum}.fastq"
    echo -e "${header}\n${seq}\n${plus}\n${qual}" > "$TMP_R_FASTQ"

    OUTPREFIX_KMC="${KMC_COUNTS_DIR}/${readname}_k${K}"
    "$KMC_DIR/kmc" -b -cs10000 -k${K} "$TMP_R_FASTQ" "$OUTPREFIX_KMC" "$KMC_COUNTS_DIR" > /dev/null 2>&1
    "$KMC_DIR/kmc_tools" -t8 transform "$OUTPREFIX_KMC" dump "${OUTPREFIX_KMC}.txt" > /dev/null 2>&1

    result=$("$SCRIPT_DIR/check_complexity" "${OUTPREFIX_KMC}.txt" "${OUTPREFIX_KMC}_complexity.log")
    if [[ "$result" == "complex" ]]; then
        if [[ -n "$header_extra" ]]; then
            echo ">$readname $header_extra" >> "$COMPLEX_FASTA"
        else
            echo ">$readname" >> "$COMPLEX_FASTA"
        fi
        echo "$seq" >> "$COMPLEX_FASTA"
    fi

    # Cleanup temporary files
    rm -f "$TMP_R_FASTQ" "${OUTPREFIX_KMC}.kmc_pre" "${OUTPREFIX_KMC}.kmc_suf"
done < "$FASTQ_FILE"

# ----------------------------
# Step 4: Run EquiRep (working directory)
# ----------------------------
info "[4/9] Running EquiRep..."
"$EQUIREP_DIR/EquiRep" "$COMPLEX_FASTA" "$WORKING_DIR/er_${INPUT_NAME}"
"$SCRIPT_DIR/update_unit_lengths.sh" "$WORKING_DIR/er_${INPUT_NAME}.fasta"

# ----------------------------
# Step 5: Run find_A (intermediate files in working dir)
# ----------------------------
info "[5/9] Running find_A..."
"$SCRIPT_DIR/find_A" \
    "$WORKING_DIR/er_${INPUT_NAME}.fasta" \
    "$WORKING_DIR/${INPUT_NAME}_predicted_A1.fasta" \
    "$WORKING_DIR/${INPUT_NAME}_predicted_Abar1.fasta" \
    "$WORKING_DIR/${INPUT_NAME}_predicted_A2.fasta" \
    "$WORKING_DIR/${INPUT_NAME}_predicted_Abar2.fasta" \
    "$WORKING_DIR/${INPUT_NAME}_predicted_AAbar.fasta" \
    "$WORKING_DIR/er_${INPUT_NAME}_rearranged.fasta" \
    "$WORKING_DIR/${INPUT_NAME}_prediction.log"

# ----------------------------
# Step 6: Generate two-copy sequences (working dir)
# ----------------------------
info "[6/9] Generating two-copy sequences..."
"$SCRIPT_DIR/make_two_copy" "$WORKING_DIR/${INPUT_NAME}_predicted_A1.fasta" "$WORKING_DIR/${INPUT_NAME}_predicted_A1_twocopy.fasta"
"$SCRIPT_DIR/make_two_copy" "$WORKING_DIR/${INPUT_NAME}_predicted_A2.fasta" "$WORKING_DIR/${INPUT_NAME}_predicted_A2_twocopy.fasta"

# ----------------------------
# Step 7: Align using minimap2 (working dir)
# ----------------------------
info "[7/9] Aligning two-copy sequences with minimap2..."
"$SCRIPT_DIR/align_minimap2.sh" "$WORKING_DIR/${INPUT_NAME}_predicted_A1_twocopy.fasta" "$WORKING_DIR/${INPUT_NAME}_predicted_A1_twocopy_align.sam" "$MINIMAP2_DIR" &>/dev/null
"$SCRIPT_DIR/align_minimap2.sh" "$WORKING_DIR/${INPUT_NAME}_predicted_A2_twocopy.fasta" "$WORKING_DIR/${INPUT_NAME}_predicted_A2_twocopy_align.sam" "$MINIMAP2_DIR" &>/dev/null

# ----------------------------
# Step 8: Calculate CIGAR (working dir)
# ----------------------------
info "[8/9] Calculating CIGAR statistics..."
"$SCRIPT_DIR/calculate_cigar" "$WORKING_DIR/${INPUT_NAME}_predicted_A1_twocopy.fasta" "$WORKING_DIR/${INPUT_NAME}_predicted_A1_twocopy_align.sam" "$WORKING_DIR/${INPUT_NAME}_predicted_A1_cigar_evaluation.log"
"$SCRIPT_DIR/calculate_cigar" "$WORKING_DIR/${INPUT_NAME}_predicted_A2_twocopy.fasta" "$WORKING_DIR/${INPUT_NAME}_predicted_A2_twocopy_align.sam" "$WORKING_DIR/${INPUT_NAME}_predicted_A2_cigar_evaluation.log"

# ----------------------------
# Step 9: Combine A1 and A2 to final A.fasta (final output prefix)
# ----------------------------
info "[9/9] Combining predictions to generate final A.fasta..."
"$SCRIPT_DIR/combine" \
    "$WORKING_DIR/${INPUT_NAME}_predicted_A1_cigar_evaluation.log" \
    "$WORKING_DIR/${INPUT_NAME}_predicted_A2_cigar_evaluation.log" \
    "$WORKING_DIR/${INPUT_NAME}_predicted_A1.fasta" \
    "$WORKING_DIR/${INPUT_NAME}_predicted_A2.fasta" \
    "${OUTPREFIX}_circRNA_seqs.fasta" "${OUTPREFIX}_circRNA_bsjs.tsv" &>/dev/null

# ----------------------------
# Cleanup working directory
# ----------------------------
info "Cleaning up temporary files..."
rm -rf "$WORKING_DIR"

info "✅ Pipeline completed successfully!"
info "Final outputs are prefixed by: $OUTPREFIX"
