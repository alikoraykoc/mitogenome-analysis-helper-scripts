#!/usr/bin/env python3
import os
import argparse
import logging
import re
from Bio import SeqIO
from Bio.SeqFeature import FeatureLocation, SeqFeature

# --- CONFIGURATION ---
# Use lowercase pairs for matching
WHITELIST_PAIRS = {("atp8", "atp6"), ("nd4", "nd4l"), ("nad4", "nad4l"), ("nad4l", "nad4")}
LENGTH_RETENTION_THRESHOLD = 0.95

# Overlap policy for NON-whitelisted pairs:
MAX_NONWHITELIST_OVERLAP_BP = 30          # if overlap > 30 bp => severe
MAX_NONWHITELIST_OVERLAP_FRAC = 0.25      # or if overlap > 25% of CDS length => severe

START_CODONS = {"ATG", "ATT", "ATA", "ATC", "GTG", "TTG"}
STOP_CODONS = {"TAA", "TAG"}

KNOWN_GENES = [
    "cox1", "cox2", "cox3",
    "atp6", "atp8",
    "nad1", "nad2", "nad3", "nad4", "nad4l", "nad5", "nad6",
    "nd1", "nd2", "nd3", "nd4", "nd4l", "nd5", "nd6",
    "cob", "cytb"
]

def setup_logging(output_dir: str):
    os.makedirs(output_dir, exist_ok=True)
    log_path = os.path.join(output_dir, "mitogenome_refinement.log")
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(message)s",
        handlers=[logging.FileHandler(log_path), logging.StreamHandler()],
    )

def get_feat_text(feat) -> str:
    q = feat.qualifiers
    terms = (
        q.get("label", [])
        + q.get("ID", [])
        + q.get("gene", [])
        + q.get("product", [])
        + q.get("note", [])
    )
    return " ".join(str(x) for x in terms).lower()

def get_gene_label(feat) -> str:
    """
    Robust gene label extractor: tries gene/label/product/note/ID and returns a normalized gene name.
    Prefers known mitochondrial PCG names (nd*, cox*, atp*, cob/cytb).
    """
    text = get_feat_text(feat)

    # Normalize common synonyms
    text = text.replace("cytochrome b", "cytb")
    text = text.replace("cytochrome-b", "cytb")
    text = text.replace("cob", "cytb")  # optional; comment out if you prefer COB label
    text = text.replace("nad", "nad")   # keep as is; we will normalize later

    # Prefer longer matches first (nad4l before nad4 etc.)
    for g in sorted(KNOWN_GENES, key=len, reverse=True):
        # accept both nad and nd variants
        variants = {g}
        if g.startswith("nad"):
            variants.add("nd" + g[3:])
        if g.startswith("nd"):
            variants.add("nad" + g[2:])

        for v in variants:
            if re.search(rf"\b{re.escape(v)}\b", text):
                # return canonical lowercase
                return v.lower()

    # Fallback: try first token-ish
    toks = re.findall(r"[A-Za-z0-9_]+", text)
    return toks[0].lower() if toks else "unknown"

def overlap_bp(loc_a, loc_b) -> int:
    a0, a1 = int(loc_a.start), int(loc_a.end)
    b0, b1 = int(loc_b.start), int(loc_b.end)
    return max(0, min(a1, b1) - max(a0, b0))

def check_orf(seq_obj, table_id=5, allow_incomplete=False) -> bool:
    """
    Table-5 ORF validity:
    - no internal stops
    - if not incomplete, must terminate with stop
    - trim to mod3 for translation stability
    """
    if len(seq_obj) < 3:
        return False
    trim_len = len(seq_obj) - (len(seq_obj) % 3)
    trimmed = seq_obj[:trim_len]
    if len(trimmed) < 3:
        return False

    aa = trimmed.translate(table=table_id)

    if "*" in aa[:-1]:
        return False
    if not allow_incomplete and (len(aa) == 0 or aa[-1] != "*"):
        return False
    return True

def find_best_orf_in_range(record, start_limit, end_limit, strand, table_id):
    """
    Find the longest valid ORF within [start_limit, end_limit], respecting strand.
    Allows incomplete stop only if ORF hits the wall exactly (candidate_end == region_len).
    """
    region_seq = record.seq[start_limit:end_limit]
    if strand == -1:
        region_seq = region_seq.reverse_complement()

    best_s, best_e, best_inc, max_len = None, None, False, 0
    region_len = len(region_seq)

    # scan all possible start positions
    for j in range(0, region_len - 2):
        codon = str(region_seq[j:j+3]).upper()
        if codon not in START_CODONS:
            continue

        curr = j + 3
        found_stop = False
        c_end = None

        while curr + 3 <= region_len:
            next_codon = str(region_seq[curr:curr+3]).upper()
            if next_codon in STOP_CODONS:
                found_stop = True
                c_end = curr + 3
                break
            curr += 3

        c_inc = False
        if not found_stop:
            # incomplete stop only if it abuts wall
            rem = (region_len - j) % 3
            if rem == 1 and str(region_seq[-1]).upper() == "T":
                c_end, c_inc, found_stop = region_len, True, True
            elif rem == 2 and str(region_seq[-2:]).upper() == "TA":
                c_end, c_inc, found_stop = region_len, True, True

        if found_stop and c_end and (c_end - j) > max_len:
            cand_seq = region_seq[j:c_end]
            if check_orf(cand_seq, table_id=table_id, allow_incomplete=c_inc):
                max_len = c_end - j
                best_inc = c_inc

                if strand == 1:
                    best_s, best_e = start_limit + j, start_limit + c_end
                else:
                    # map back to genomic coords
                    best_s, best_e = end_limit - c_end, end_limit - j

    return best_s, best_e, best_inc, max_len

def process_batch(input_dir, output_dir, table_id=5):
    setup_logging(output_dir)
    files = [f for f in os.listdir(input_dir) if f.endswith((".gb", ".gbk"))]

    if not files:
        logging.warning("No .gb/.gbk files found in input directory.")
        return

    for filename in files:
        sample_id = os.path.splitext(filename)[0]
        logging.info(f"--- Processing: {sample_id} ---")

        record = SeqIO.read(os.path.join(input_dir, filename), "genbank")
        seq_len = len(record.seq)

        # 1) CATEGORIZE FEATURES
        trnas, rrnas, pcgs = [], [], []
        for f in record.features:
            txt = get_feat_text(f)
            if f.type == "tRNA":
                trnas.append(f)
            elif f.type == "CDS":
                pcgs.append(f)
            elif f.type == "rRNA" or any(x in txt for x in ["rrn", "16s", "12s"]):
                f.type = "rRNA"
                rrnas.append(f)

        trnas.sort(key=lambda x: int(x.location.start))
        pcgs.sort(key=lambda x: int(x.location.start))

        # 2) APPLY rRNA & D-LOOP RULES
        rrnL = next((r for r in rrnas if any(x in get_feat_text(r) for x in ["rrnl", "16s"])), None)
        rrnS = next((r for r in rrnas if any(x in get_feat_text(r) for x in ["rrns", "12s"])), None)
        trnV = next((t for t in trnas if any(x in get_feat_text(t) for x in ["trnv", "val"])), None)

        if rrnL:
            leu_trnas = [t for t in trnas if any(x in get_feat_text(t) for x in ["leu", "trnl"])]
            if leu_trnas and trnV:
                trnL1 = min(leu_trnas, key=lambda t: abs(int(t.location.end) - int(rrnL.location.start)))
                rrnL.location = FeatureLocation(trnL1.location.end, trnV.location.start, strand=-1)
                logging.info(f"  [rRNA] Snapped rrnL between {get_feat_text(trnL1)} and {get_feat_text(trnV)}")

        if trnV and rrnS:
            rrnS.location = FeatureLocation(trnV.location.end, rrnS.location.end, strand=-1)
            logging.info(f"  [rRNA] Snapped rrnS after {get_feat_text(trnV)}")

        new_features = list(trnas) + list(rrnas)

        if rrnS:
            d_loop = SeqFeature(FeatureLocation(rrnS.location.end, seq_len, strand=1), type="misc_feature")
            d_loop.qualifiers["label"] = ["Putative_D_loop_Control_Region"]
            new_features.append(d_loop)

        # Helper: quick access to nearest tRNA walls around a coordinate
        def upstream_trna_end_before(pos: int) -> int:
            return max([int(t.location.end) for t in trnas if int(t.location.end) <= pos] + [0])

        def downstream_trna_start_after(pos: int) -> int:
            return min([int(t.location.start) for t in trnas if int(t.location.start) >= pos] + [seq_len])

        # 3) CDS REFINEMENT
        last_fixed_end = 0

        for i, cds in enumerate(pcgs):
            label = get_gene_label(cds)  # robust label

            orig_start = int(cds.location.start)
            orig_end = int(cds.location.end)
            orig_len = int(cds.location.end) - int(cds.location.start)
            strand = cds.location.strand

            prev_label = get_gene_label(pcgs[i-1]) if i > 0 else None
            next_label = get_gene_label(pcgs[i+1]) if i < len(pcgs) - 1 else None

            # whitelist check (order matters)
            is_exception = False
            if prev_label and (prev_label, label) in WHITELIST_PAIRS:
                is_exception = True
            if next_label and (label, next_label) in WHITELIST_PAIRS:
                is_exception = True

            # Determine deterministic walls for this CDS (0-overlap policy)
            up_wall = max(last_fixed_end, upstream_trna_end_before(orig_start))
            down_wall = min(downstream_trna_start_after(orig_end),
                            int(pcgs[i+1].location.start) if i < len(pcgs) - 1 else seq_len)

            # Safety: if walls collapse or invert, skip ORF search and flag
            if down_wall <= up_wall:
                cds.qualifiers.setdefault("note", []).append("WALLS_INVALID: down_wall<=up_wall; manual review")
                logging.warning(f"  [{label.upper()}] INVALID WALLS (up={up_wall}, down={down_wall}). Keeping original.")
                new_features.append(cds)
                last_fixed_end = int(cds.location.end)
                continue

            # If exception, keep as-is (biological overlap allowed for those pairs)
            if is_exception:
                new_features.append(cds)
                last_fixed_end = int(cds.location.end)
                logging.info(f"  [{label.upper()}] Whitelisted pair overlap allowed. Kept original.")
                continue

            # Compute overlap severity with the previous placed CDS (if any)
            prev_cds = next((f for f in reversed(new_features) if f.type == "CDS"), None)
            ov = overlap_bp(cds.location, prev_cds.location) if prev_cds else 0
            ov_frac = ov / max(1, orig_len)

            # If current CDS violates walls or overlaps previous CDS, attempt fix
            violates = (orig_start < up_wall) or (orig_end > down_wall) or (ov > 0)

            if violates:
                best_s, best_e, is_inc, new_len = find_best_orf_in_range(record, up_wall, down_wall, strand, table_id)

                if best_s is not None and new_len >= orig_len * LENGTH_RETENTION_THRESHOLD:
                    cds.location = FeatureLocation(best_s, best_e, strand=strand)
                    logging.info(f"  [{label.upper()}] Fixed within walls. New length={new_len} (orig={orig_len}).")
                else:
                    # If severe overlap (non-whitelist), DO NOT silently accept
                    if (ov > MAX_NONWHITELIST_OVERLAP_BP) or (ov_frac > MAX_NONWHITELIST_OVERLAP_FRAC):
                        cds.qualifiers.setdefault("note", []).append(
                            f"SEVERE_OVERLAP_NONWHITELIST: ov={ov}bp ({ov_frac:.1%}); manual review"
                        )
                        logging.warning(
                            f"  [{label.upper()}] SEVERE overlap non-whitelist (ov={ov}bp, {ov_frac:.1%}). FLAGGED."
                        )
                    else:
                        logging.info(f"  [{label.upper()}] Minor overlap or wall violation could not be fixed; preserved.")
            else:
                logging.info(f"  [{label.upper()}] No overlap/wall violation detected.")

            new_features.append(cds)
            last_fixed_end = int(cds.location.end)

        # Final sort and write
        record.features = sorted(new_features, key=lambda x: int(x.location.start))
        out_path = os.path.join(output_dir, f"{sample_id}_refined.gb")
        SeqIO.write(record, out_path, "genbank")
        logging.info(f"Saved: {out_path}")

def main():
    parser = argparse.ArgumentParser(description="Batch refine mitogenome GenBank annotations.")
    parser.add_argument("-i", "--input", required=True, help="Input directory with .gb/.gbk files")
    parser.add_argument("-o", "--output", required=True, help="Output directory for refined files")
    parser.add_argument("--table", type=int, default=5, help="NCBI translation table (default 5)")
    args = parser.parse_args()
    process_batch(args.input, args.output, table_id=args.table)

if __name__ == "__main__":
    main()
