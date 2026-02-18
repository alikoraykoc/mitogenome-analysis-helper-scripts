#!/usr/bin/env python3

import copy
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Seq import Seq

def check_orf(seq_obj, table=5, allow_incomplete=False):
    if len(seq_obj) < 3:
        return False
        
    trim_len = len(seq_obj) - (len(seq_obj) % 3)
    trimmed = seq_obj[:trim_len]
    
    if len(trimmed) < 3:
        return False
    
    translation = trimmed.translate(table=table)
    
    if "*" in translation[:-1]:
        return False
        
    if not allow_incomplete and (len(translation) == 0 or translation[-1] != "*"):
        return False
        
    return True

def get_feat_text(feat):
    quals = feat.qualifiers
    terms = quals.get("label", []) + quals.get("ID", []) + quals.get("gene", []) + quals.get("product", []) + quals.get("standard_name", [])
    return " ".join(terms).lower()

def process_mitogenome(input_gb, output_gb):
    record = SeqIO.read(input_gb, "genbank")
    seq_len = len(record.seq)
    
    trnas, pcgs, rrnas = [], [], []

    for feat in record.features:
        text = get_feat_text(feat)
        
        if feat.type == "tRNA":
            trnas.append(feat)
        elif feat.type == "CDS": 
            pcgs.append(feat)
        elif feat.type == "rRNA" or "rrn" in text or "16s" in text or "12s" in text:
            feat.type = "rRNA"
            rrnas.append(feat)
            
    trnas.sort(key=lambda f: int(f.location.start))
    new_features = list(trnas) + list(rrnas)

    print("\n" + "="*65)
    print(" MITOGENOME CURATION & QC REPORT (STRICT ZERO-OVERLAP)")
    print("="*65)

    # 1. Akıllı rRNA Çapaları (İsim yerine mesafeye göre trnL1 bulma)
    rrnL = next((r for r in rrnas if "rrnl" in get_feat_text(r) or "16s" in get_feat_text(r) or "l-rrna" in get_feat_text(r)), None)
    rrnS = next((r for r in rrnas if "rrns" in get_feat_text(r) or "12s" in get_feat_text(r) or "s-rrna" in get_feat_text(r)), None)
    trnV = next((t for t in trnas if "trnv" in get_feat_text(t) or "val" in get_feat_text(t)), None)
    
    trnL1 = None
    if rrnL:
        leu_trnas = [t for t in trnas if "leu" in get_feat_text(t) or "trnl" in get_feat_text(t)]
        if leu_trnas:
            # rrnL'ye fiziksel olarak en yakın olan tRNA-Leu'yu trnL1 kabul et
            trnL1 = min(leu_trnas, key=lambda t: abs(int(t.location.end) - int(rrnL.location.start)))

    if trnL1 and trnV and rrnL:
        rrnL.location = FeatureLocation(trnL1.location.end, trnV.location.start, strand=-1)
        print("[rRNA] rrnL (16S) boundaries strictly snapped to trnL1 and trnV.")
    else:
        print("[rRNA] WARNING: Anchors for rrnL missing. Preserved original boundaries.")

    if trnV and rrnS:
        rrnS.location = FeatureLocation(trnV.location.end, rrnS.location.end, strand=-1)
        print("[rRNA] rrnS (12S) start boundary strictly snapped to trnV.")
    else:
        print("[rRNA] WARNING: trnV anchor missing. Preserved original rrnS boundaries.")
        
    if rrnS:
        d_loop = SeqFeature(FeatureLocation(rrnS.location.end, seq_len, strand=1), type="misc_feature")
        d_loop.qualifiers["label"] = ["Putative_D_loop_Control_Region"]
        new_features.append(d_loop)
        print("[D-LOOP] Created spanning from rrnS terminus to sequence end.")
    
    print("-" * 65)

    pcg_labels = [p.qualifiers.get("label", p.qualifiers.get("gene", ["Unknown"]))[0].lower() for p in pcgs]
    if len(pcgs) != 13:
        print(f"[CRITICAL QC WARNING] Expected 13 CDS, but identified {len(pcgs)}.")
    
    pcgs.sort(key=lambda f: int(f.location.start))
    last_fixed_pcg_end = 0 
    
    # 2. PCG Duvarları (Midpoint Mantığı ile Kusursuz Çakışma Engelleme)
    for i, pcg in enumerate(pcgs):
        orig_start = int(pcg.location.start)
        orig_end = int(pcg.location.end)
        strand = pcg.location.strand
        quals = pcg.qualifiers
        label = quals.get("label", quals.get("gene", quals.get("product", ["Unknown_CDS"])))[0]
        
        # [BUG FIX]: Duvarları orig_end/start yerine genin orta noktasına (midpoint) göre belirle!
        midpoint = (orig_start + orig_end) / 2
        
        upstream_trna_end = max([int(t.location.end) for t in trnas if int(t.location.end) < midpoint] + [0])
        upstream_pcg_end = last_fixed_pcg_end 
        min_pos = max(upstream_trna_end, upstream_pcg_end)
        
        downstream_trna_start = min([int(t.location.start) for t in trnas if int(t.location.start) > midpoint] + [seq_len])
        downstream_pcg_start = int(pcgs[i+1].location.start) if i < len(pcgs)-1 else seq_len
        max_pos = min(downstream_trna_start, downstream_pcg_start)

        wall_type = "tRNA" if max_pos == downstream_trna_start else "PCG"

        region_seq = record.seq[min_pos:max_pos]
        if strand == -1:
            region_seq = region_seq.reverse_complement()
            
        region_len = len(region_seq)
        search_start = 0
        search_end = region_len - 2

        best_start = None
        best_end = None
        is_incomplete = False
        longest_valid_orf_len = 0 

        for j in range(search_start, search_end):
            codon = str(region_seq[j:j+3]).upper()
            
            if codon in ["ATG", "ATT", "ATA", "ATC", "GTG", "TTG"]:
                current_idx = j + 3
                found_stop = False
                candidate_end = None
                
                while current_idx + 3 <= region_len:
                    next_codon = str(region_seq[current_idx:current_idx+3]).upper()
                    if next_codon in ["TAA", "TAG"]:
                        found_stop = True
                        candidate_end = current_idx + 3
                        break
                    current_idx += 3
                
                candidate_incomplete = False
                
                if not found_stop:
                    remainder = (region_len - j) % 3
                    if remainder == 1 and str(region_seq[-1]).upper() == "T":
                        candidate_end = region_len
                        candidate_incomplete = True
                        found_stop = True
                    elif remainder == 2 and str(region_seq[-2:]).upper() == "TA":
                        candidate_end = region_len
                        candidate_incomplete = True
                        found_stop = True

                if found_stop and candidate_end:
                    candidate_seq = region_seq[j:candidate_end]
                    
                    if check_orf(candidate_seq, table=5, allow_incomplete=candidate_incomplete):
                        orf_length = candidate_end - j
                        
                        if orf_length > longest_valid_orf_len:
                            longest_valid_orf_len = orf_length
                            
                            if strand == 1:
                                best_start = min_pos + j
                                best_end = min_pos + candidate_end
                            else:
                                best_start = max_pos - candidate_end
                                best_end = max_pos - j
                                
                            is_incomplete = candidate_incomplete

        if best_start is not None and best_end is not None:
            pcg.location = FeatureLocation(best_start, best_end, strand=strand)
            last_fixed_pcg_end = int(pcg.location.end)
            
            seq_len_qc = abs(best_end - best_start)
            mod_qc = "OK" if seq_len_qc % 3 == 0 or is_incomplete else "FAIL"
            inc_qc = f"(Incomplete Stop: T/TA abutting {wall_type})" if is_incomplete else ""
            
            print(f"[{label.upper():<6}] Length: {seq_len_qc:<4}bp | Frame: {mod_qc} {inc_qc:<40} | Status: FIXED")
            
            if is_incomplete:
                pcg.qualifiers["note"] = [f"Incomplete stop codon abutting downstream {wall_type} (0-bp overlap enforced)."]
        else:
            print(f"[{label.upper():<6}] WARNING: No valid ORF identified without overlap! Manual curation required.")
            
            fallback_start = max(orig_start, min_pos)
            fallback_end = min(orig_end, max_pos)
            pcg.location = FeatureLocation(fallback_start, fallback_end, strand=strand)
            last_fixed_pcg_end = int(pcg.location.end)
            
            pcg.qualifiers["note"] = ["WARNING: Could not resolve strict zero-overlap constraint. Putative biological overlap or frameshift requires manual review."]

        new_features.append(pcg)

    new_features.sort(key=lambda f: int(f.location.start))
    record.features = new_features
    SeqIO.write(record, output_gb, "genbank")
    print("-" * 65)
    print(f"[SUCCESS] Curated annotation saved to: {output_gb}")

if __name__ == "__main__":
    INPUT_FILE = "file/path.gb" 
    OUTPUT_FILE = "file/path.gb"
    
    process_mitogenome(INPUT_FILE, OUTPUT_FILE)
