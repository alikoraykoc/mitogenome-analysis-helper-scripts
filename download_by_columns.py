#!/usr/bin/env python3
"""
Download GenBank sequences column-by-column from an Excel sheet and
merge each column into a single FASTA named after the gene (column name).

Behavior requested by user:
- Read the Excel **column by column**.
- Start reading columns **after** the "species" column.
- For each target column (a gene), download all accessions in that column.
- Write them into **one FASTA file per column**, named "{gene}.fasta".
- Each sequence header should be **genus_species** where genus/species
  are parsed from the "species" column of the same row.
- After finishing a column, print an English message like "Finished column: {gene}".

Usage:
  python download_by_columns.py --excel "Tetrigidae genbank.xlsx" --sheet "Sheet1" --outdir out \
      --email you@example.com [--api_key YOUR_NCBI_API_KEY]

Notes:
  - Requires: pandas, biopython
  - Respects simple rate-limiting for NCBI (sleep ~0.35s/request).
"""

import argparse
import os
import re
import sys
import time
from typing import List, Optional
import pandas as pd
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord
from io import StringIO

SAFE_FILENAME_CHARS = r"[^A-Za-z0-9._+-]"

def sanitize_filename(name: str) -> str:
    return re.sub(SAFE_FILENAME_CHARS, "_", name.strip())

def parse_species_cell(cell: str) -> Optional[tuple]:
    """
    From a 'species' cell like 'Tettigidea lateralis (Say, 1824)'
    return ('tettigidea', 'lateralis').
    If parsing fails, return None.
    """
    if not isinstance(cell, str):
        return None
    # Take the first two words that look like a binomial
    tokens = cell.strip().split()
    if len(tokens) < 2:
        return None
    genus = re.sub(r"[^A-Za-z]", "", tokens[0])
    species = re.sub(r"[^A-Za-z]", "", tokens[1])
    if not genus or not species:
        return None
    return (genus.lower(), species.lower())

def split_accessions(cell) -> List[str]:
    """
    Split a cell into possible accession tokens.
    Accept comma/semicolon/whitespace separated accessions.
    """
    if cell is None:
        return []
    if isinstance(cell, float):
        # NaN
        return []
    if not isinstance(cell, str):
        cell = str(cell)
    # normalize separators
    cell = cell.strip()
    if not cell:
        return []
    # Replace common separators with spaces
    cell = re.sub(r"[;,]", " ", cell)
    parts = [p.strip() for p in cell.split() if p.strip()]
    # filter plausible accession-like tokens (allow letters, digits, underscore, dot)
    accs = [p for p in parts if re.fullmatch(r"[A-Za-z0-9_.]+", p)]
    return accs

def fetch_fasta_by_accession(acc: str, email: str, api_key: Optional[str] = None, db: str = "nuccore") -> SeqRecord:
    """
    Fetch a nucleotide sequence in FASTA by accession.
    Returns a SeqRecord.
    """
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    # efetch FASTA
    with Entrez.efetch(db=db, id=acc, rettype="fasta", retmode="text") as handle:
        fasta_txt = handle.read().strip()
    # Parse to SeqRecord
    records = list(SeqIO.parse(StringIO(fasta_txt), "fasta"))
    if not records:
        raise ValueError(f"No FASTA returned for {acc}")
    return records[0]

def main():
    ap = argparse.ArgumentParser(description="Download GenBank sequences column-by-column from Excel.")
    ap.add_argument("--excel", required=True, help="Path to Excel file (.xlsx)")
    ap.add_argument("--sheet", default=None, help="Sheet name (default: first sheet)")
    ap.add_argument("--outdir", required=True, help="Output directory for per-gene FASTA files")
    ap.add_argument("--email", required=True, help="Email for NCBI Entrez")
    ap.add_argument("--api_key", default=None, help="NCBI API key (optional)")
    ap.add_argument("--start_after_header", default="species",
                    help='Start after the last column whose header matches this (case-insensitive substring). Default: "species"')
    ap.add_argument("--sleep", type=float, default=0.35, help="Sleep seconds between NCBI requests (default 0.35)")
    ap.add_argument("--db", default="nuccore", help="NCBI database (default nuccore)")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    # Read Excel
    df = pd.read_excel(args.excel, sheet_name=args.sheet)
    # Identify 'species' column
    species_col_candidates = [c for c in df.columns if str(c).strip().lower() == "species"]
    if not species_col_candidates:
        print("ERROR: Could not find a 'species' column in the Excel.", file=sys.stderr)
        sys.exit(1)
    species_col = species_col_candidates[0]

    # Determine start column index: after header containing args.start_after_header
    start_after = args.start_after_header.strip().lower()
    start_idx = -1
    for i, col in enumerate(df.columns):
        if start_after in str(col).strip().lower():
            start_idx = i
    if start_idx == -1:
        print(f"WARNING: No column header containing '{args.start_after_header}' was found. Starting from the first column.", file=sys.stderr)
        start_idx = -1

    # Target columns are those strictly after start_idx
    target_cols = [c for c in df.columns[start_idx+1:]]

    # Also skip obvious non-gene descriptor columns if they happen to come after
    skip_names = {"family", "subfamily", "tribe", "genus"}
    filtered = []
    for c in target_cols:
        if str(c).strip().lower() in skip_names:
            continue
        filtered.append(c)
    target_cols = filtered

    if not target_cols:
        print("No target gene columns detected after the specified header.", file=sys.stderr)
        sys.exit(1)

    print(f"Detected target columns (genes): {', '.join(str(c) for c in target_cols)}")

    for col in target_cols:
        gene = str(col).strip()
        gene_safe = sanitize_filename(gene)
        out_fasta = os.path.join(args.outdir, f"{gene_safe}.fasta")
        n_written = 0

        with open(out_fasta, "w") as fout:
            # iterate rows
            for idx, cell in df[col].items():
                accs = split_accessions(cell)
                if not accs:
                    continue

                # derive genus/species from species column of same row
                sp = df.at[idx, species_col] if idx in df.index else None
                parsed = parse_species_cell(sp)
                if not parsed:
                    # fallback if species column is not parsable
                    genus_species = f"unknown_{idx}"
                else:
                    genus, species = parsed
                    genus_species = f"{genus}_{species}"

                for acc in accs:
                    try:
                        rec = fetch_fasta_by_accession(acc, email=args.email, api_key=args.api_key, db=args.db)
                        # build custom header
                        header = f"{genus_species}"
                        # create new SeqRecord with our header
                        new_rec = SeqRecord(rec.seq, id=header, description="")
                        SeqIO.write(new_rec, fout, "fasta")
                        n_written += 1
                        time.sleep(args.sleep)
                    except Exception as e:
                        print(f"WARNING: failed to fetch {acc} for row {idx}, gene {gene}: {e}", file=sys.stderr)
                        # continue to next accession

        print(f"Finished column: {gene} | wrote {n_written} sequences -> {out_fasta}")

if __name__ == "__main__":
    main()
