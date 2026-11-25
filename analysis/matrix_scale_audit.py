"""Summarise the actual processed inputs without exporting case-level values."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from core_final import load_apgi, load_cptac, load_fudan, load_ruijin


def main(root, output):
    sources = [
        ("Fudan", lambda: load_fudan(root / "13045_2022_1384_MOESM25_ESM.xlsx"), "log2(x+1) intensity"),
        ("CPTAC", lambda: load_cptac(root / "cptac_orig/proteomics_gene_level_MD_abundance_tumor.cct",
                                     root / "cptac_orig/proteomics_gene_level_MD_abundance_normal.cct"),
         "released common-reference TMT log2 abundance"),
        ("APGI", lambda: load_apgi(root / "apgi/protein_matrix.csv", root / "hgnc_complete_set.txt"),
         "released DIA/MaxLFQ analysed scale"),
        ("Ruijin", lambda: load_ruijin(root / "node_ruijin/20221227_PDAC_PRO_exp.xlsx"),
         "released processed abundance; upstream denominator and imputation unspecified"),
    ]
    summaries = []
    for label, load, scale in sources:
        result = load()
        raw = result["abundance_X"]
        finite = np.isfinite(raw)
        zero = finite & (raw == 0)
        values = raw[finite]
        summaries.append(dict(cohort=label, scale=scale, release_protein_rows=raw.shape[0],
                              release_samples=raw.shape[1], finite=int(finite.sum()), missing=int((~finite).sum()),
                              finite_exact_zero=int(zero.sum()), finite_negative=int((values < 0).sum()),
                              min=values.min(), median=np.median(values), max=values.max(),
                              labels="before complete-pair and nomination filters; CPTAC uses intersecting protein rows"))
        print(label, summaries[-1], flush=True)
    pd.DataFrame(summaries).to_csv(output, sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    main(args.inputs, args.output)
