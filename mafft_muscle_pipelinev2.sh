#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# MAFFT & MUSCLE ALIGNMENT PIPELINE
# Description: Adjusts sequence direction with MAFFT, removes gaps, and 
#              generates final alignments using MUSCLE.
# ==============================================================================

# --- CONFIGURATION ---
# Define input directory containing raw FASTA files
INDIR="/Users/koray/Desktop/phylogeny_mitogenome/nov-2025/run7_onlymtg/dataset"

# Define output directories for each pipeline stage
MAFFT_DIR="01_mafft"
UNGAP_DIR="02_mafft_ungapped"
MUSCLE_DIR="03_muscle"

# Define executable names
MAFFT_BIN="mafft"
MUSCLE_BIN="muscle"  # Uses MUSCLE5 syntax (-align / -output)

# --- INITIALIZATION ---
mkdir -p "$MAFFT_DIR" "$UNGAP_DIR" "$MUSCLE_DIR"

echo "[INFO] Initializing Alignment Pipeline"
echo "[INFO] Input Directory  : $INDIR"
echo "[INFO] MAFFT Output     : $MAFFT_DIR"
echo "[INFO] Ungapped Output  : $UNGAP_DIR"
echo "[INFO] MUSCLE Output    : $MUSCLE_DIR"
echo "--------------------------------------------------------------------------------"

# --- PROCESSING LOOP ---
shopt -s nullglob
for fasta in "$INDIR"/*.fasta; do
    base=$(basename "$fasta")
    stem="${base%.fasta}"

    echo "[INFO] Processing: $base"

    # Step 1: Direction adjustment and initial alignment via MAFFT
    mafft_out="$MAFFT_DIR/${stem}_mafft.fasta"
    echo "       -> Running MAFFT (Direction Adjustment)..."
    "$MAFFT_BIN" --auto --adjustdirection "$fasta" > "$mafft_out" 2>/dev/null

    # Step 2: Remove gaps and MAFFT's '_R_' prefix from reverse-complemented sequences
    ungap_out="$UNGAP_DIR/${stem}_mafft_ungapped.fasta"
    echo "       -> Formatting sequences (Removing gaps and '_R_' prefixes)..."
    awk 'BEGIN{OFS=""} /^>/{sub(">_R_", ">"); print; next} {gsub("-","",$0); print}' \
        "$mafft_out" > "$ungap_out"

    # Step 3: Final multiple sequence alignment via MUSCLE
    muscle_out="$MUSCLE_DIR/${stem}_muscle.fasta"
    echo "       -> Running MUSCLE..."
    "$MUSCLE_BIN" -align "$ungap_out" -output "$muscle_out" >/dev/null 2>&1

    echo "       -> Completed: $muscle_out"
    echo "--------------------------------------------------------------------------------"
done

echo "[INFO] Pipeline execution finished successfully."