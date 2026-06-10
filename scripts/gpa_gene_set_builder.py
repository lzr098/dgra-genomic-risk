#!/usr/bin/env python3
"""
gene_set_builder.py - 动态构建表型相关基因集

功能: 根据患者表型关键词，通过 OMIM 本地数据库和 HPO 数据构建候选基因集。
      可作为大模型推理后的验证/补充工具，也可独立使用。

用法:
    python3 gene_set_builder.py \
        --phenotypes "joubert,polydactyly,epilepsy" \
        --omim-db ~/.workbuddy/data/omim/omim.db \
        --output genes.txt \
        --max-genes 50

输出:
    每行一个基因符号的文本文件
"""

import argparse
import sqlite3
import sys
import re


# HPO ID to gene mapping (simplified — can be expanded or linked to HPO API)
# This is a minimal built-in mapping for common phenotypes.
# For production use, download full HPO annotations from:
# https://hpo.jax.org/app/download/annotation
HPO_GENE_MAP = {
    # Ciliopathy / Joubert
    "joubert": ["OFD1", "CEP290", "NPHP1", "AHI1", "CC2D2A", "RPGRIP1L",
                "TMEM67", "TMEM216", "MKS1", "B9D1", "B9D2", "ARL13B",
                "INPP5E", "TCTN1", "TCTN2", "TCTN3", "CSPP1", "KIF7",
                "KIAA0556", "KIAA0753", "PDE6D", "NPHP3", "NPHP4",
                "SDCCAG8", "WDR35", "WDR19", "WDR60", "IFT80", "IFT140",
                "IFT172", "TTC21B", "C5orf42", "ARMC9", "CASC8"],

    # Polydactyly
    "polydactyly": ["GLI3", "SHH", "LMBR1", "MBTS1", "ZRS", "IQCE",
                    "MUTYH", "OFD1", "KIAA0753"],

    # Epilepsy
    "epilepsy": ["SCN1A", "SCN2A", "SCN8A", "KCNQ2", "KCNQ3", "STXBP1",
                 "CDKL5", "PCDH19", "ARX", "SLC2A1", "SPTAN1", "TSC1",
                 "TSC2", "MECP2", "FOXG1", "TCF4", "EHMT1", "KCNB1",
                 "GABRA1", "GABRB3", "GRIN2A", "GRIN2B", "CACNA1A",
                 "CHD2", "DYRK1A", "SYNGAP1", "SCN1B", "SCN9A",
                 "KCNT1", "KCNA2", "KCNMA1", "SLC6A1", "HCN1",
                 "GABRG2", "DEPDC5", "NPRL2", "NPRL3", "PRICKLE1",
                 "PRICKLE2", "LGI1", "CASK", "PNKP", "POLG",
                 "ALDH7A1", "PNPO", "PLPBP"],

    # Intellectual disability / developmental delay
    "intellectual_disability": ["MECP2", "CDKL5", "FOXG1", "TCF4", "EHMT1",
                                 "ARID1B", "KDM6B", "KDM6A", "MEF2C", "SATB2",
                                 "DYRK1A", "SYNGAP1", "KCNB1", "CHD2", "SETD5",
                                 "POGZ", "ANKRD11", "KAT6A", "KMT2A", "KMT2D",
                                 "KMT2E", "SMARCA2", "SMARCA4", "SMARCB1",
                                 "SMARCE1", "ARID1A", "ARID2", "PBRM1",
                                 "BRWD3", "HUWE1", "UPF3B", "ZDHHC9",
                                 "IL1RAPL1", "TSPAN7", "OPHN1", "PAK3",
                                 "SLC16A2", "ZIC2", "ZIC3", "NRXN1",
                                 "CNTNAP2", "SHANK2", "SHANK3", "NLGN3",
                                 "NLGN4X", "NRXN2", "NRXN3", "GRIN2A",
                                 "GRIN2B", "GRIA3", "CACNA1C", "CACNA1D",
                                 "CACNA1H", "SCN2A", "SCN8A", "KCNQ2",
                                 "KCNQ3", "KCNB1", "KCNT1", "HCN1"],

    # Craniofacial / micrognathia
    "craniofacial": ["FGFR1", "FGFR2", "FGFR3", "TWIST1", "EFNB1", "MSX1",
                     "MSX2", "PAX3", "PAX6", "PAX9", "ALX1", "ALX3",
                     "ALX4", "PRRX1", "PRRX2", "HOXA2", "HOXA1", "HOXB1",
                     "SIX1", "SIX2", "SIX3", "SIX5", "EYA1", "EYA4",
                     "TCOF1", "POLR1C", "POLR1D", "SF3B4", "KMT2D",
                     "KDM6A", "ARID1B", "SATB2", "GLI2", "GLI3",
                     "SHH", "PTCH1", "SUFU", "OFD1"],

    # Brain malformation / neuronal migration
    "brain_malformation": ["TUBB2B", "TUBA1A", "TUBB3", "TUBG1", "DYNC1H1",
                           "KIF5C", "KIF2A", "LIS1", "DCX", "ARX", "RELN",
                           "VLDLR", "WDR62", "CDK5RAP2", "CEP152", "STIL",
                           "CENPJ", "ASPM", "MCPH1", "WDR35", "WDR19",
                           "OFD1", "CEP290", "NPHP1", "AHI1", "CC2D2A",
                           "TMEM67", "FOXC1", "ZIC1", "ZIC2", "ZIC3",
                           "ZIC4", "ZIC5", "PAX6", "EMX2", "SOX2"],

    # Enlarged cisterna magna / Dandy-Walker
    "cisterna_magna": ["OFD1", "CEP290", "NPHP1", "AHI1", "CC2D2A",
                       "TMEM67", "TMEM216", "MKS1", "B9D1", "B9D2",
                       "FOXC1", "ZIC1", "ZIC4", "WDR35", "WDR19"],

    # Retinitis pigmentosa (ciliopathy-related)
    "retinitis_pigmentosa": ["OFD1", "CEP290", "NPHP1", "NPHP3", "NPHP4",
                             "NPHP5", "NPHP6", "RPGR", "RPGRIP1", "RPGRIP1L",
                             "USH2A", "ABCA4", "RHO", "PDE6A", "PDE6B",
                             "CNGB1", "CNGA1", "MERTK", "PRPH2", "RPE65",
                             "CRB1", "EYS", "CERKL", "SNRNP200", "SEMA4A"],

    # Polycystic kidney (ciliopathy-related)
    "polycystic_kidney": ["PKD1", "PKD2", "PKHD1", "OFD1", "NPHP1", "NPHP3",
                          "NPHP4", "NPHP5", "NPHP6", "NPHP7", "NPHP8",
                          "NPHP9", "NPHP10", "NPHP11", "NPHP12", "NPHP13",
                          "NPHP14", "NPHP15", "NPHP16", "NPHP17", "NPHP18",
                          "NPHP19", "NPHP20", "HNF1B", "UMOD", "REN",
                          "ACE", "AGT", "AGTR1", "DYNC2H1", "WDR35",
                          "WDR19", "IFT80", "IFT140", "IFT172", "TTC21B"],

    # Simpson-Golabi-Behmel / overgrowth
    "simpson_golabi_behmel": ["GPC3", "GPC4", "OFD1"],

    # Ciliopathy (generic)
    "ciliopathy": ["OFD1", "CEP290", "NPHP1", "NPHP3", "NPHP4", "NPHP5",
                   "NPHP6", "NPHP7", "NPHP8", "NPHP9", "NPHP10", "NPHP11",
                   "NPHP12", "NPHP13", "NPHP14", "NPHP15", "NPHP16",
                   "NPHP17", "NPHP18", "NPHP19", "NPHP20", "AHI1",
                   "CC2D2A", "RPGRIP1L", "TMEM67", "TMEM216", "MKS1",
                   "B9D1", "B9D2", "ARL13B", "INPP5E", "TCTN1",
                   "TCTN2", "TCTN3", "CSPP1", "KIF7", "KIAA0556",
                   "KIAA0753", "PDE6D", "SDCCAG8", "WDR35", "WDR19",
                   "WDR60", "IFT80", "IFT140", "IFT172", "TTC21B",
                   "C5orf42", "ARMC9", "CASC8", "DYNC2H1"],

    # Developmental delay (generic)
    "developmental_delay": ["MECP2", "CDKL5", "FOXG1", "TCF4", "EHMT1",
                            "ARID1B", "KDM6B", "KDM6A", "MEF2C", "SATB2",
                            "DYRK1A", "SYNGAP1", "KCNB1", "CHD2", "SETD5",
                            "POGZ", "ANKRD11", "KAT6A", "KMT2D", "KMT2A",
                            "SMARCA2", "SMARCA4", "SMARCB1", "HUWE1",
                            "NRXN1", "SHANK2", "SHANK3", "NLGN3", "NLGN4X",
                            "GRIN2A", "GRIN2B", "CACNA1C", "SCN2A",
                            "KCNQ2", "KCNQ3", "HNRNPU", "STAG1", "STAG2",
                            "NIPBL", "RAD21", "SMC1A", "SMC3", "HDAC8",
                            "TBR1", "TRIO", "PPP2R5D", "CSNK2A1"],
}


def search_omim(omim_db_path, keywords):
    """从 OMIM 本地数据库搜索相关基因"""
    genes = set()
    if not omim_db_path:
        return genes

    try:
        conn = sqlite3.connect(omim_db_path)
        cursor = conn.cursor()

        # Try different table schemas
        tables = []
        cursor.execute("SELECT name FROM sqlite_master WHERE type='table'")
        for row in cursor.fetchall():
            tables.append(row[0])

        # Search by title and text fields
        for keyword in keywords:
            keyword_lower = keyword.lower()

            # Try 'omim' table (most common schema)
            if 'omim' in tables:
                try:
                    cursor.execute(
                        "SELECT gene_symbol, symbols FROM omim WHERE LOWER(title) LIKE ? OR LOWER(text) LIKE ?",
                        (f'%{keyword_lower}%', f'%{keyword_lower}%')
                    )
                    for row in cursor.fetchall():
                        if row[0]:
                            genes.add(row[0])
                        # Parse secondary symbols
                        if row[1]:
                            for sym in str(row[1]).split(','):
                                genes.add(sym.strip())
                except Exception:
                    pass

            # Try 'genes' table
            if 'genes' in tables:
                try:
                    cursor.execute(
                        "SELECT gene_symbol FROM genes WHERE LOWER(gene_name) LIKE ?",
                        (f'%{keyword_lower}%',)
                    )
                    for row in cursor.fetchall():
                        if row[0]:
                            genes.add(row[0])
                except Exception:
                    pass

        conn.close()
    except Exception as e:
        print(f"WARNING: OMIM query failed: {e}", file=sys.stderr)

    return genes


def build_gene_set(phenotypes, omim_db=None, max_genes=100):
    """
    构建表型相关基因集。

    策略：
    1. 从内置 HPO 映射获取基因
    2. 从 OMIM 本地数据库补充
    3. 去重、排序、截断
    """
    all_genes = set()

    # 1. Built-in HPO mapping
    for pheno in phenotypes:
        pheno_clean = pheno.lower().strip().replace(' ', '_').replace('-', '_')
        if pheno_clean in HPO_GENE_MAP:
            all_genes.update(HPO_GENE_MAP[pheno_clean])
            print(f"  [{pheno_clean}] built-in: {len(HPO_GENE_MAP[pheno_clean])} genes")
        else:
            # Try partial match
            matched = False
            for key, genes in HPO_GENE_MAP.items():
                if pheno_clean in key or key in pheno_clean:
                    all_genes.update(genes)
                    print(f"  [{pheno_clean}] -> matched '{key}': {len(genes)} genes")
                    matched = True
                    break
            if not matched:
                print(f"  [{pheno_clean}] no built-in mapping found")

    # 2. OMIM supplement
    if omim_db:
        omim_genes = search_omim(omim_db, phenotypes)
        if omim_genes:
            all_genes.update(omim_genes)
            print(f"  [OMIM] added {len(omim_genes)} genes")

    # Sort and limit
    gene_list = sorted(all_genes)
    if len(gene_list) > max_genes:
        print(f"WARNING: Gene set truncated from {len(gene_list)} to {max_genes}", file=sys.stderr)
        gene_list = gene_list[:max_genes]

    return gene_list


def main():
    parser = argparse.ArgumentParser(description='Build phenotype-driven gene set')
    parser.add_argument('--phenotypes', required=True,
                        help='Comma-separated phenotype keywords (e.g., "joubert,polydactyly,epilepsy")')
    parser.add_argument('--omim-db',
                        help='Path to OMIM SQLite database (optional)')
    parser.add_argument('--output', required=True, help='Output gene list file')
    parser.add_argument('--max-genes', type=int, default=100,
                        help='Maximum number of genes to return')

    args = parser.parse_args()

    phenotypes = [p.strip() for p in args.phenotypes.split(',')]
    print(f"Building gene set for phenotypes: {phenotypes}")

    genes = build_gene_set(phenotypes, omim_db=args.omim_db, max_genes=args.max_genes)

    print(f"\nTotal unique genes: {len(genes)}")

    with open(args.output, 'w') as f:
        for gene in genes:
            f.write(gene + '\n')

    print(f"Gene list written to: {args.output}")
    print(f"First 20 genes: {', '.join(genes[:20])}")


if __name__ == '__main__':
    main()
