#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# FINAL HYBRID ALIGNMENT PIPELINE: MACSE (PCGs) & MUSCLE (RNAs)
# Description: Adjusts sequence direction, cleans sequence labels, aligns PCGs 
#              via codon-aware MACSE (Table 5), and aligns RNAs via MUSCLE.
#              Replaces MACSE '!' flags with '-' for downstream compatibility.
# ==============================================================================

# --- CONFIGURATION ---
INDIR=""

MAFFT_DIR="01_mafft"
UNGAP_DIR="02_mafft_ungapped"
ALIGN_DIR="03_aligned"

MAFFT_BIN="mafft"
MUSCLE_BIN="muscle"
MACSE_BIN="macse"

# --- INITIALIZATION ---
mkdir -p "$MAFFT_DIR" "$UNGAP_DIR" "$ALIGN_DIR"

echo "[INFO] Starting Final Hybrid Alignment Pipeline..."
echo "[INFO] Input Directory : $INDIR"
echo "[INFO] Final Outputs   : $ALIGN_DIR"
echo "--------------------------------------------------------------------------------"

# Regex pattern to identify Protein-Coding Genes (PCGs) - case-insensitive
PCG_PATTERN="cox|atp|nad|nd[1-6]|cytb|cob"

# --- PROCESSING LOOP ---
shopt -s nullglob
for fasta in "$INDIR"/*.fasta; do
    base=$(basename "$fasta")
    stem="${base%.fasta}"

    echo "[INFO] Processing: $base"

    # Step 1: Adjust sequence direction using MAFFT
    mafft_out="$MAFFT_DIR/${stem}_mafft.fasta"
    "$MAFFT_BIN" --auto --adjustdirection "$fasta" > "$mafft_out" 2>/dev/null

    # Step 2: Clean alignment gaps ('-') and reverse-complement prefixes ('_R_') using AWK
    ungap_out="$UNGAP_DIR/${stem}_mafft_ungapped.fasta"
    awk 'BEGIN{OFS=""} /^>/{sub(">_R_", ">"); print; next} {gsub("-","",$0); print}' \
        "$mafft_out" > "$ungap_out"

    # Step 3: Determine alignment strategy based on gene type
    align_out="$ALIGN_DIR/${stem}_aligned.fasta"

    if echo "$stem" | grep -qiE "$PCG_PATTERN"; then
        echo "       -> Running MACSE (Codon-based alignment for PCG)..."
        
        # Align using MACSE with the Invertebrate Mitochondrial genetic code (Translation Table 5)
        "$MACSE_BIN" -prog alignSequences -gc_def 5 -seq "$ungap_out" -out_NT "$align_out" >/dev/null 2>&1
        
        # Replace MACSE-specific '!' frameshift flags with standard '-' gap characters for IQ-TREE
        sed 's/!/-/g' "$align_out" > "${align_out}.tmp" && mv "${align_out}.tmp" "$align_out"
        
    else
        echo "       -> Running MUSCLE (Standard nucleotide alignment for RNA)..."
        "$MUSCLE_BIN" -align "$ungap_out" -output "$align_out" >/dev/null 2>&1
    fi

    echo "       -> Completed: $align_out"
    echo "--------------------------------------------------------------------------------"
done

echo "[SUCCESS] All genes aligned and cleaned for downstream phylogenetic analysis."
