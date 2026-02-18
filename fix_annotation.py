import os
import csv
import argparse
import logging
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import FeatureLocation, SeqFeature

# --- GLOBAL CONFIGURATION ---
WHITELIST_PAIRS = [("atp8", "atp6"), ("nd4", "nd4l"), ("nad4", "nad4l")]

# Standard Invertebrate Mitochondrial Code (Table 5) Mapping for RSCU
AMINO_ACID_MAP = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'M', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': 'W', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def setup_logging(output_dir):
    log_file = os.path.join(output_dir, "pipeline_log.txt")
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s [%(levelname)s] %(message)s',
        handlers=[logging.FileHandler(log_file), logging.StreamHandler()]
    )

def calculate_skews(sequence):
    s = sequence.upper()
    a, t, g, c = s.count("A"), s.count("T"), s.count("G"), s.count("C")
    at_skew = (a - t) / (a + t) if (a + t) > 0 else 0
    gc_skew = (g - c) / (g + c) if (g + c) > 0 else 0
    at_content = (a + t) / len(s) if len(s) > 0 else 0
    return round(at_content, 4), round(at_skew, 4), round(gc_skew, 4)

def calculate_rscu(pcg_sequences):
    combined_seq = "".join(pcg_sequences).upper()
    codons = [combined_seq[i:i+3] for i in range(0, len(combined_seq)-len(combined_seq)%3, 3)]
    counts = Counter(codons)
    
    aa_counts = {}
    for codon, count in counts.items():
        if codon in AMINO_ACID_MAP:
            aa = AMINO_ACID_MAP[codon]
            if aa not in aa_counts: aa_counts[aa] = {}
            aa_counts[aa][codon] = count

    rscu_data = []
    for aa, synonymous_codons in aa_counts.items():
        total_aa_count = sum(synonymous_codons.values())
        n_i = len([c for c, a in AMINO_ACID_MAP.items() if a == aa])
        for codon, count in synonymous_codons.items():
            rscu = count / (total_aa_count / n_i)
            rscu_data.append({"AminoAcid": aa, "Codon": codon, "Count": count, "RSCU": round(rscu, 4)})
    return rscu_data

def process_batch(input_dir, output_dir, table_id):
    setup_logging(output_dir)
    dirs = {
        "gb": os.path.join(output_dir, "Fixed_GenBanks"),
        "nt": os.path.join(output_dir, "Extracted_Genes_NT"),
        "aa": os.path.join(output_dir, "Extracted_Genes_AA")
    }
    for d in dirs.values(): os.makedirs(d, exist_ok=True)

    summary_stats = []
    all_rscu_profiles = []
    gb_files = [f for f in os.listdir(input_dir) if f.endswith((".gb", ".gbk"))]

    logging.info(f"Starting batch process for {len(gb_files)} files using Table {table_id}")

    for gb_file in gb_files:
        sample_id = os.path.splitext(gb_file)[0]
        logging.info(f"--- Processing Sample: {sample_id} ---")
        
        try:
            record = SeqIO.read(os.path.join(input_dir, gb_file), "genbank")
        except Exception as e:
            logging.error(f"Failed to read {gb_file}: {e}"); continue

        seq_len = len(record.seq)
        trnas, pcgs, rrnas = [], [], []
        for feat in record.features:
            q = feat.qualifiers
            txt = " ".join(q.get("label", []) + q.get("gene", []) + q.get("product", [])).lower()
            if feat.type == "tRNA": trnas.append(feat)
            elif feat.type == "CDS": pcgs.append(feat)
            elif "rrn" in txt or feat.type == "rRNA": feat.type = "rRNA"; rrnas.append(feat)

        pcgs.sort(key=lambda f: int(f.location.start))
        trnas.sort(key=lambda f: int(f.location.start))
        
        last_end = 0
        final_feats = list(trnas) + list(rrnas)
        sample_pcg_seqs = []

        for i, pcg in enumerate(pcgs):
            lbl = " ".join(pcg.qualifiers.get("label", pcg.qualifiers.get("gene", ["PCG"]))).split()[0].upper()
            
            # Whitelist Allowance Logic
            up_allow = 7 if (i > 0 and (pcgs[i-1].qualifiers.get("label", [""])[0].lower(), lbl.lower()) in WHITELIST_PAIRS) else 0
            down_allow = 7 if (i < len(pcgs)-1 and (lbl.lower(), pcgs[i+1].qualifiers.get("label", [""])[0].lower()) in WHITELIST_PAIRS) else 0
            
            mid = (int(pcg.location.start) + int(pcg.location.end)) / 2
            min_pos = max(max([int(t.location.end) for t in trnas if int(t.location.end) < mid] + [0]), last_end - up_allow)
            max_pos = min(min([int(t.location.start) for t in trnas if int(t.location.start) > mid] + [seq_len]), 
                         (int(pcgs[i+1].location.start) + down_allow) if i < len(pcgs)-1 else seq_len)

            best_orf_seq, best_s, best_e, is_inc = None, None, None, False
            region = record.seq[min_pos:max_pos]
            if pcg.location.strand == -1: region = region.reverse_complement()
            
            # ORF Search
            for j in range(len(region)-2):
                if str(region[j:j+3]).upper() in ["ATG", "ATT", "ATA", "ATC", "GTG", "TTG"]:
                    curr, found_stop, c_end = j+3, False, None
                    while curr+3 <= len(region):
                        if str(region[curr:curr+3]).upper() in ["TAA", "TAG"]: found_stop, c_end = True, curr+3; break
                        curr += 3
                    
                    c_inc = False
                    if not found_stop:
                        if (len(region)-j)%3 == 1 and str(region[-1]).upper() == "T": c_end, c_inc, found_stop = len(region), True, True
                        elif (len(region)-j)%3 == 2 and str(region[-2:]).upper() == "TA": c_end, c_inc, found_stop = len(region), True, True

                    if found_stop and c_end:
                        orf_cand = region[j:c_end]
                        if "*" not in orf_cand[:len(orf_cand)-(len(orf_cand)%3)].translate(table=table_id)[:-1]:
                            if pcg.location.strand == 1: best_s, best_e = min_pos + j, min_pos + c_end
                            else: best_s, best_e = max_pos - c_end, max_pos - j
                            best_orf_seq, is_inc = orf_cand, c_inc; break

            if best_s is not None:
                pcg.location = FeatureLocation(best_s, best_e, strand=pcg.location.strand)
                last_end = int(pcg.location.end)
                sample_pcg_seqs.append(str(best_orf_seq))
                
                # Stats & Fasta
                at_c, at_s, gc_s = calculate_skews(best_orf_seq)
                summary_stats.append({
                    "Sample": sample_id, "Gene": lbl, "Length": len(best_orf_seq),
                    "AT_Content": at_c, "AT_Skew": at_s, "GC_Skew": gc_s,
                    "Start": str(best_orf_seq[:3]), "Stop": str(best_orf_seq[-3:]) if not is_inc else "INC"
                })
                
                with open(os.path.join(dirs["nt"], f"{lbl}.fasta"), "a") as f:
                    f.write(f">{sample_id}\n{best_orf_seq}\n")
                with open(os.path.join(dirs["aa"], f"{lbl}.fasta"), "a") as f:
                    aa = best_orf_seq[:len(best_orf_seq)-(len(best_orf_seq)%3)].translate(table=table_id)
                    f.write(f">{sample_id}\n{aa}\n")
            else:
                logging.warning(f"No valid ORF for {lbl} in {sample_id}!")
            final_feats.append(pcg)

        # RSCU Processing
        if sample_pcg_seqs:
            rscu_profile = calculate_rscu(sample_pcg_seqs)
            for entry in rscu_profile: entry["Sample"] = sample_id; all_rscu_profiles.append(entry)

        record.features = final_feats
        SeqIO.write(record, os.path.join(dirs["gb"], f"{sample_id}_fixed.gb"), "genbank")

    # CSV Exports
    with open(os.path.join(output_dir, "Summary_Statistics.csv"), "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=summary_stats[0].keys())
        writer.writeheader(); writer.writerows(summary_stats)
    with open(os.path.join(output_dir, "RSCU_Profiles.csv"), "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=all_rscu_profiles[0].keys())
        writer.writeheader(); writer.writerows(all_rscu_profiles)

    logging.info("Batch Processing Complete. Check output directory for results.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Professional Mitogenome Curation & Extraction Pipeline")
    parser.add_argument("-i", "--input", required=True, help="Input directory containing .gb files")
    parser.add_argument("-o", "--output", required=True, help="Output directory for results")
    parser.add_argument("-t", "--table", type=int, default=5, help="Translation table (default: 5)")
    args = parser.parse_args()
    process_batch(args.input, args.output, args.table)
