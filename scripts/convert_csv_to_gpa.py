#!/usr/bin/env python3
"""
Convert annotation CSV → GPA-compatible TSV
"""
import os
import csv
import sys

def detect_encoding(filepath):
    """Try multiple encodings to find the correct one."""
    # Read raw bytes to detect
    with open(filepath, 'rb') as f:
        raw = f.read(50000)
    
    # Try encodings in order of likelihood for Chinese bioinformatics files
    encodings = ['gb18030', 'gbk', 'utf-8', 'latin1', 'cp1252']
    
    for enc in encodings:
        try:
            decoded = raw.decode(enc)
            # Check for common bioinformatics terms that should be readable
            if 'CHROM' in decoded or 'gene' in decoded.lower() or 'format_GT' in decoded:
                # Verify GT field looks correct
                lines = decoded.split('\n')
                if len(lines) > 1:
                    header = lines[0]
                    gt_idx = header.split(',').index('format_GT') if 'format_GT' in header else -1
                    if gt_idx >= 0 and len(lines) > 1:
                        first_row = lines[1].split(',')
                        if len(first_row) > gt_idx:
                            gt_val = first_row[gt_idx]
                            if gt_val in ('0/1', '1/1', './.', '0/0', '0|1', '1|1'):
                                print(f"[ENCODING] Detected {enc} (GT field valid: {gt_val})")
                                return enc
        except Exception:
            pass
    
    # Fallback: try each encoding on the full file
    for enc in encodings:
        try:
            with open(filepath, 'r', encoding=enc) as f:
                header = f.readline()
                first = f.readline()
                parts = first.split(',')
                if len(parts) > 71:  # format_GT is at index 71
                    gt_val = parts[71].strip()
                    if gt_val in ('0/1', '1/1', './.', '0/0', '0|1', '1|1'):
                        print(f"[ENCODING] Detected {enc} (GT field valid: {gt_val})")
                        return enc
        except Exception:
            pass
    
    print("[ENCODING] WARNING: Could not detect encoding, falling back to gb18030")
    return 'gb18030'


# Consequence → IMPACT mapping (VEP-style)
IMPACT_MAP = {
    'transcript_ablation': 'HIGH',
    'splice_acceptor_variant': 'HIGH',
    'splice_donor_variant': 'HIGH',
    'stop_gained': 'HIGH',
    'frameshift_variant': 'HIGH',
    'stop_lost': 'HIGH',
    'start_lost': 'HIGH',
    'transcript_amplification': 'HIGH',
    'inframe_insertion': 'MODERATE',
    'inframe_deletion': 'MODERATE',
    'missense_variant': 'MODERATE',
    'protein_altering_variant': 'MODERATE',
    'splice_region_variant': 'LOW',
    'splice_donor_region_variant': 'LOW',
    'splice_polypyrimidine_tract_variant': 'LOW',
    'incomplete_terminal_codon_variant': 'LOW',
    'start_retained_variant': 'LOW',
    'stop_retained_variant': 'LOW',
    'synonymous_variant': 'LOW',
    'coding_sequence_variant': 'MODERATE',
    'mature_miRNA_variant': 'LOW',
    '5_prime_UTR_variant': 'MODIFIER',
    '3_prime_UTR_variant': 'MODIFIER',
    'non_coding_transcript_exon_variant': 'MODIFIER',
    'intron_variant': 'MODIFIER',
    'NMD_transcript_variant': 'MODIFIER',
    'non_coding_transcript_variant': 'MODIFIER',
    'upstream_gene_variant': 'MODIFIER',
    'downstream_gene_variant': 'MODIFIER',
    'TFBS_ablation': 'MODIFIER',
    'TFBS_amplification': 'MODIFIER',
    'TF_binding_site_variant': 'MODIFIER',
    'regulatory_region_ablation': 'MODIFIER',
    'regulatory_region_amplification': 'MODIFIER',
    'feature_elongation': 'MODIFIER',
    'regulatory_region_variant': 'MODIFIER',
    'feature_truncation': 'MODIFIER',
    'intergenic_variant': 'MODIFIER',
}


def infer_impact(consequence):
    """Infer VEP IMPACT from consequence term."""
    if not consequence:
        return 'MODIFIER'
    first = consequence.split(',')[0].split(';')[0].strip()
    return IMPACT_MAP.get(first, 'MODIFIER')


def convert(input_csv, output_tsv):
    encoding = detect_encoding(input_csv)
    
    with open(input_csv, 'r', encoding=encoding) as fin, \
         open(output_tsv, 'w', encoding='utf-8', newline='') as fout:
        reader = csv.DictReader(fin)
        
        out_cols = [
            'CHROM', 'POS', 'REF', 'ALT', 'GENE',
            'Consequence', 'IMPACT', 'HGVSc', 'HGVSp',
            'GT', 'DP', 'GQ', 'VAF',
        ]
        writer = csv.DictWriter(fout, fieldnames=out_cols, delimiter='\t')
        writer.writeheader()
        
        count = 0
        gt_errors = 0
        for row in reader:
            gt = row.get('format_GT', './.').strip()
            vaf = ''
            
            # Validate GT field
            if gt in ('0/1', '1/1', '0/0', './.', '0|1', '1|1', '0|0'):
                if gt == '1/1' or gt == '1|1':
                    vaf = '1.0'
                elif gt == '0/1' or gt == '0|1':
                    vaf = '0.5'
                elif gt == '0/0' or gt == '0|0':
                    vaf = '0.0'
            else:
                gt_errors += 1
                if gt_errors <= 5:
                    print(f"[GT WARNING] Invalid GT value: '{gt}' (gene={row.get('gene','')}, pos={row.get('pos','')})")
                gt = './.'
            
            consequence = row.get('consequence', '')
            impact = infer_impact(consequence)
            
            out_row = {
                'CHROM': row.get('chrom', ''),
                'POS': row.get('pos', ''),
                'REF': row.get('ref', ''),
                'ALT': row.get('alt', ''),
                'GENE': row.get('gene', row.get('vep_default_gene', '')),
                'Consequence': consequence,
                'IMPACT': impact,
                'HGVSc': row.get('hgvsc', ''),
                'HGVSp': row.get('hgvsp', ''),
                'GT': gt,
                'DP': row.get('format_DP', ''),
                'GQ': row.get('format_GQ', ''),
                'VAF': vaf,
            }
            writer.writerow(out_row)
            count += 1
        
        print(f"Converted {count} variants → {output_tsv}")
        if gt_errors > 0:
            print(f"[WARNING] {gt_errors} variants had invalid GT values (set to ./.)")


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python3 convert_csv_to_gpa.py <input.csv> <output.tsv>")
        sys.exit(1)
    convert(sys.argv[1], sys.argv[2])
