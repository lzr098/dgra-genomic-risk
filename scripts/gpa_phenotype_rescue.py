#!/usr/bin/env python3
"""
phenotype_rescue.py - 表型驱动VCF候选变异救援搜索 (通用版)

功能：在VEP注释后的VCF中，根据动态构建的表型相关基因集搜索候选变异。
      基因集不是硬编码的，而是通过外部数据库或用户输入动态确定。

用法：
    # 方式1: 直接提供基因列表
    python3 phenotype_rescue.py \
        --vcf P002.vep.vcf.gz \
        --genes "OFD1,CEP290,WDR35,GLI3" \
        --output P002.rescue.tsv \
        --patient-sex male

    # 方式2: 从文件读取基因列表
    python3 phenotype_rescue.py \
        --vcf P002.vep.vcf.gz \
        --gene-list genes.txt \
        --output P002.rescue.tsv

输出列：
    Priority, Gene, Chrom, Pos, Ref, Alt, GT, DP, GQ, gnomAD_AF,
    Impact, Consequence, HGVSc, HGVSp, Exon, ClinVar, Feature,
    Rescue_Reason
"""

import argparse
import gzip
import re
import sys


def parse_vep_csq_format(vcf_header_lines):
    """从VCF header解析CSQ字段格式"""
    for line in vcf_header_lines:
        if "CSQ=Format:" in line:
            match = re.search(r'Format: ([^"]+)', line)
            if match:
                return match.group(1).split('|')
    # 默认字段列表（VEP标准）
    return [
        "Allele", "Consequence", "IMPACT", "SYMBOL", "Gene", "Feature_type",
        "Feature", "BIOTYPE", "EXON", "INTRON", "HGVSc", "HGVSp",
        "cDNA_position", "CDS_position", "Protein_position", "Amino_acids",
        "Codons", "Existing_variation", "DISTANCE", "STRAND", "FLAGS",
        "SYMBOL_SOURCE", "HGNC_ID", "gnomADe_AF", "gnomADe_AFR_AF",
        "gnomADe_AMR_AF", "gnomADe_ASJ_AF", "gnomADe_EAS_AF",
        "gnomADe_FIN_AF", "gnomADe_MID_AF", "gnomADe_NFE_AF",
        "gnomADe_REMAINING_AF", "gnomADe_SAS_AF", "CLIN_SIG",
        "SOMATIC", "PHENO"
    ]


def determine_priority(variant, patient_sex):
    """确定救援优先级"""
    reasons = []
    priority = 3  # 默认低优先级

    gt = variant.get('gt', '')
    af = variant.get('gnomad_af')
    chrom = variant.get('chrom', '')
    impact = variant.get('impact', '')

    is_x_linked = chrom in ('chrX', 'X')
    is_male = patient_sex == 'male'

    # P0: X-linked + male + any non-ref GT = hemizygous (most important!)
    if is_x_linked and is_male and gt in ('1/1', '1|1', '0/1', '0|1', '1|0'):
        priority = min(priority, 1)
        reasons.append("X-linked_hemizygous_male")

    # P0: 纯合/半合子 + 罕见
    elif gt in ('1/1', '1|1'):
        if is_x_linked and is_male:
            priority = min(priority, 1)
            reasons.append("X-linked_hemizygous_male")
        else:
            priority = min(priority, 1)
            reasons.append("homozygous")

    # P1: 极低频率
    if af is not None and af != '':
        try:
            af_val = float(af)
            if af_val < 0.0001:
                priority = min(priority, 1)
                reasons.append(f"ultra_rare_AF={af_val}")
            elif af_val < 0.001:
                priority = min(priority, 2)
                reasons.append(f"very_rare_AF={af_val}")
        except:
            pass
    elif af == '' or af is None:
        priority = min(priority, 2)
        reasons.append("not_in_gnomAD")

    # HIGH impact 升级
    if impact == 'HIGH':
        priority = min(priority, 1)
        reasons.append("HIGH_impact")

    return priority, ";".join(reasons) if reasons else "rescue_search"


def search_vcf(vcf_path, target_genes, min_impact="MODERATE", max_af=0.01,
               patient_sex=None):
    """在VCF中搜索目标基因的候选变异"""

    impact_levels = {"HIGH": 3, "MODERATE": 2, "LOW": 1, "MODIFIER": 0}
    min_level = impact_levels.get(min_impact, 2)

    header_lines = []
    csq_fields = None
    csq_idx = {}
    results = []

    target_genes = set(g.strip() for g in target_genes)

    opener = gzip.open if str(vcf_path).endswith('.gz') else open

    with opener(vcf_path, 'rt') as f:
        for line in f:
            if line.startswith('##'):
                header_lines.append(line)
                continue
            if line.startswith('#CHROM'):
                header_lines.append(line)
                csq_fields = parse_vep_csq_format(header_lines)
                csq_idx = {f: i for i, f in enumerate(csq_fields)}
                break

        if not csq_fields:
            print("ERROR: Could not parse VEP CSQ format from header", file=sys.stderr)
            return []

        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            chrom, pos, id_, ref, alt = parts[0], parts[1], parts[2], parts[3], parts[4]
            info = parts[7]
            fmt = parts[8] if len(parts) > 8 else ''
            sample = parts[9] if len(parts) > 9 else ''

            csq_match = re.search(r'CSQ=([^;]+)', info)
            if not csq_match:
                continue

            csq_entries = csq_match.group(1).split(',')
            for csq_str in csq_entries:
                csq_parts = csq_str.split('|')
                if len(csq_parts) < len(csq_fields):
                    continue

                gene = csq_parts[csq_idx.get('SYMBOL', 3)]
                if gene not in target_genes:
                    continue

                impact = csq_parts[csq_idx.get('IMPACT', 2)]
                if impact_levels.get(impact, 0) < min_level:
                    continue

                af = csq_parts[csq_idx.get('gnomADe_AF', 23)]
                if af and af != '':
                    try:
                        if float(af) >= max_af:
                            continue
                    except:
                        pass

                gt = dp = gq = ''
                if fmt and sample:
                    for k, v in zip(fmt.split(':'), sample.split(':')):
                        if k == 'GT':
                            gt = v
                        elif k == 'DP':
                            dp = v
                        elif k == 'GQ':
                            gq = v

                variant = {
                    'gene': gene,
                    'chrom': chrom,
                    'pos': pos,
                    'ref': ref,
                    'alt': alt,
                    'gt': gt,
                    'dp': dp,
                    'gq': gq,
                    'gnomad_af': af if af else None,
                    'impact': impact,
                    'consequence': csq_parts[csq_idx.get('Consequence', 1)],
                    'hgvsc': csq_parts[csq_idx.get('HGVSc', 10)],
                    'hgvsp': csq_parts[csq_idx.get('HGVSp', 11)],
                    'exon': csq_parts[csq_idx.get('EXON', 8)],
                    'clinvar': csq_parts[csq_idx.get('CLIN_SIG', 32)],
                    'feature': csq_parts[csq_idx.get('Feature', 6)],
                }

                priority, reason = determine_priority(variant, patient_sex)
                variant['priority'] = priority
                variant['reason'] = reason
                results.append(variant)
                break

    return results


def main():
    parser = argparse.ArgumentParser(description='Phenotype-driven VCF rescue search (generic)')
    parser.add_argument('--vcf', required=True, help='VEP-annotated VCF (.vcf or .vcf.gz)')
    parser.add_argument('--genes', help='Comma-separated gene symbols')
    parser.add_argument('--gene-list', help='File with one gene per line')
    parser.add_argument('--output', required=True, help='Output TSV file')
    parser.add_argument('--min-impact', default='MODERATE',
                        choices=['HIGH', 'MODERATE', 'LOW'],
                        help='Minimum VEP IMPACT to include')
    parser.add_argument('--max-af', type=float, default=0.01,
                        help='Maximum gnomAD allele frequency')
    parser.add_argument('--patient-sex', choices=['male', 'female'],
                        help='Patient sex (affects X-linked hemizygosity scoring)')

    args = parser.parse_args()

    # Build gene list
    target_genes = set()
    if args.genes:
        target_genes.update(g.strip() for g in args.genes.split(',') if g.strip())
    if args.gene_list:
        with open(args.gene_list) as f:
            target_genes.update(line.strip() for line in f if line.strip())

    if not target_genes:
        print("ERROR: No target genes provided. Use --genes or --gene-list", file=sys.stderr)
        sys.exit(1)

    print(f"Target genes: {len(target_genes)}")
    print(f"Searching VCF: {args.vcf}")

    results = search_vcf(
        args.vcf,
        target_genes,
        min_impact=args.min_impact,
        max_af=args.max_af,
        patient_sex=args.patient_sex
    )

    print(f"Found {len(results)} candidate variants")

    with open(args.output, 'w') as out:
        header = ['Priority', 'Gene', 'Chrom', 'Pos', 'Ref', 'Alt', 'GT', 'DP', 'GQ',
                  'gnomAD_AF', 'Impact', 'Consequence', 'HGVSc', 'HGVSp', 'Exon',
                  'ClinVar', 'Feature', 'Rescue_Reason']
        out.write('\t'.join(header) + '\n')

        for v in sorted(results, key=lambda x: (x['priority'], x['gene'])):
            row = [
                str(v['priority']),
                v['gene'], v['chrom'], v['pos'], v['ref'], v['alt'],
                v['gt'], v['dp'], v['gq'],
                str(v['gnomad_af']) if v['gnomad_af'] else '',
                v['impact'], v['consequence'], v['hgvsc'], v['hgvsp'], v['exon'],
                v['clinvar'], v['feature'], v['reason']
            ]
            out.write('\t'.join(row) + '\n')

    print(f"Output written to: {args.output}")

    p1 = sum(1 for r in results if r['priority'] == 1)
    p2 = sum(1 for r in results if r['priority'] == 2)
    p3 = sum(1 for r in results if r['priority'] == 3)
    print(f"Priority breakdown: P1={p1}, P2={p2}, P3={p3}")


if __name__ == '__main__':
    main()
