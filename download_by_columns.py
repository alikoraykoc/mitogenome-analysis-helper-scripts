#!/usr/bin/env python3

"""
Excel to FASTA GenBank Downloader

This script reads a formatted Excel spreadsheet containing species names and GenBank 
accession numbers, and downloads the corresponding sequences from NCBI. It merges 
each column into a single FASTA file named after the gene, with sequence headers 
automatically formatted as 'genus_species'.

REQUIREMENTS:
    pip install pandas biopython openpyxl requests

EXPECTED EXCEL FORMAT:
    - Must contain a column named "species" (case-insensitive).
    - Cells in the "species" column should contain the binomial name (e.g., "Tettigidea lateralis").
    - All columns to the right of the "species" column are treated as target genes.
    - Cells under the gene columns should contain the GenBank accession number(s).
    - Can safely handle empty cells and non-accession text notes.

BASIC USAGE:
    python download_genbank.py --excel "Tetrigidae genbank.xlsx" --outdir output_fastas --email your@email.com

ADVANCED USAGE (Specific Gene & GFF Annotations):
    python download_genbank.py --excel "Tetrigidae genbank.xlsx" --outdir output_fastas --email your@email.com --column COX1 --gff --log download_cox1.log

ARGUMENTS:
    --excel     : Path to the input Excel (.xlsx) file (Required)
    --outdir    : Directory to save the output FASTA and GFF files (Required)
    --email     : Your NCBI email address for Entrez access (Required)
    --api_key   : Your NCBI API key to increase rate limits (Optional)
    --sheet     : Specific sheet name to read (Optional, defaults to first sheet)
    --column    : Only process this specific column/gene (Optional)
    --gff       : Also download GFF3 annotation files alongside FASTAs (Optional)
    --log       : Path to the log file (Defaults to 'download.log')
    --sleep     : Seconds to pause between NCBI requests (Defaults to 0.35s)

NOTES:
    - Automatically retries failed NCBI requests up to 3 times if the connection drops.
    - Skips downloading GFF files if they already exist in the output folder, making it easy to resume stopped jobs.
    - Pauses between downloads to respect NCBI rate limits and keep your IP from getting blocked.
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
from datetime import datetime
import requests

SAFE_FILENAME_CHARS = r"[^A-Za-z0-9._+-]"

# ================= LOGGING =================

def log_message(logfile, msg):
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    with open(logfile, "a") as f:
        f.write(f"[{ts}] {msg}\n")


# ================= UTILITIES =================

def sanitize_filename(name: str) -> str:
    return re.sub(SAFE_FILENAME_CHARS, "_", name.strip())

def parse_species_cell(cell: str) -> Optional[tuple]:
    if not isinstance(cell, str):
        return None
    tokens = cell.strip().split()
    if len(tokens) < 2:
        return None
    genus = re.sub(r"[^A-Za-z]", "", tokens[0])
    species = re.sub(r"[^A-Za-z]", "", tokens[1])
    if not genus or not species:
        return None
    return (genus.lower(), species.lower())

def split_accessions(cell) -> List[str]:
    if cell is None:
        return []
    if isinstance(cell, float):
        return []
    if not isinstance(cell, str):
        cell = str(cell)
    cell = re.sub(r"[;,]", " ", cell.strip())
    parts = [p.strip() for p in cell.split() if p.strip()]
    return [p for p in parts if re.fullmatch(r"[A-Za-z0-9_.]+", p)]


# ================= FETCHERS =================

def fetch_fasta_by_accession(acc: str, email: str, api_key=None, db="nuccore"):
    Entrez.email = email
    if api_key:
        Entrez.api_key = api_key

    for attempt in range(3):
        try:
            with Entrez.efetch(db=db, id=acc, rettype="fasta", retmode="text") as handle:
                txt = handle.read().strip()

            recs = list(SeqIO.parse(StringIO(txt), "fasta"))
            if not recs:
                raise ValueError(f"No FASTA for {acc}")
            return recs[0]
        except Exception as e:
            if attempt == 2:
                raise e
            time.sleep(2)


def fetch_gff(acc: str, outdir: str, logfile: str, email: str, api_key: Optional[str] = None, genus_species: str = "unknown"):
    outfile = os.path.join(outdir, f"{genus_species}.gff3")

    if os.path.exists(outfile):
        log_message(logfile, f"SKIP GFF exists: {acc}")
        return

    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "nuccore", 
        "id": acc, 
        "rettype": "gff3", 
        "retmode": "text",
        "email": email
    }
    if api_key:
        params["api_key"] = api_key

    for attempt in range(3):
        try:
            r = requests.get(url, params=params, timeout=15)
            if r.status_code == 200 and r.text.strip():
                with open(outfile, "w") as f:
                    f.write(r.text)
                log_message(logfile, f"GFF downloaded: {acc}")
                return
        except requests.exceptions.RequestException:
            pass
        time.sleep(2)
        
    log_message(logfile, f"GFF FAILED: {acc} after 3 attempts")


# ================= MAIN =================

def main():
    ap = argparse.ArgumentParser(description="Download GenBank sequences from Excel.")
    ap.add_argument("--excel", required=True, help="Path to Excel file")
    ap.add_argument("--sheet", default=None, help="Sheet name")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--email", required=True, help="Your NCBI email")
    ap.add_argument("--api_key", default=None, help="Your NCBI API key")

    ap.add_argument("--column", help="Process only this column")
    ap.add_argument("--gff", action="store_true", help="Also download GFF3")
    ap.add_argument("--log", default="download.log", help="Path to log file")

    ap.add_argument("--sleep", type=float, default=0.35, help="Sleep time between requests")
    args = ap.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    log_message(args.log, "=== DOWNLOAD START ===")

    # Read Excel
    df = pd.read_excel(args.excel, sheet_name=args.sheet)

    # Identify species column
    species_cols = [c for c in df.columns if str(c).strip().lower() == "species"]
    if not species_cols:
        print("ERROR: Could not find a 'species' column in the Excel.", file=sys.stderr)
        sys.exit(1)
    species_col = species_cols[0]

    # COLUMN SELECTION
    if args.column:
        if args.column not in df.columns:
            print(f"ERROR: Column '{args.column}' not found.", file=sys.stderr)
            sys.exit(1)
        target_cols = [args.column]
    else:
        start_idx = list(df.columns).index(species_col)
        target_cols = df.columns[start_idx+1:]

    # Filter out obvious non-gene columns
    skip_names = {"family", "subfamily", "tribe", "genus"}
    target_cols = [c for c in target_cols if str(c).strip().lower() not in skip_names]

    if not target_cols:
        print("No target gene columns detected.", file=sys.stderr)
        sys.exit(1)

    print("Processing columns:", list(target_cols))

    for col in target_cols:
        gene = str(col).strip()
        out_fasta = os.path.join(args.outdir, f"{sanitize_filename(gene)}.fasta")
        written = 0

        with open(out_fasta, "w") as fout:
            for idx, cell in df[col].items():

                accs = split_accessions(cell)
                if not accs:
                    continue

                sp = df.at[idx, species_col] if idx in df.index else None
                parsed = parse_species_cell(sp)
                header = f"{parsed[0]}_{parsed[1]}" if parsed else f"unknown_{idx}"

                for acc in accs:
                    try:
                        rec = fetch_fasta_by_accession(acc, args.email, args.api_key)
                        new_rec = SeqRecord(rec.seq, id=header, description="")
                        SeqIO.write(new_rec, fout, "fasta")
                        written += 1
                        log_message(args.log, f"FASTA OK: {acc} for {gene}")

                        if args.gff:
                            fetch_gff(acc, args.outdir, args.log, args.email, args.api_key, header)

                    except Exception as e:
                        log_message(args.log, f"ERROR {acc}: {e}")
                    
                    finally:
                        # Sleep runs regardless of success or failure to prevent API bans
                        time.sleep(args.sleep)

        print(f"Finished column: {gene} | wrote {written} sequences")
        log_message(args.log, f"Finished column: {gene} | wrote {written} sequences")

    log_message(args.log, "=== DOWNLOAD END ===")

if __name__ == "__main__":
    main()