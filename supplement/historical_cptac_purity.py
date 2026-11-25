"""Tissue-composition analysis in CPTAC-PDAC.

Q: is part of the protein-level support explained by tumour purity / stroma content?
- validate protein-based epithelial/stromal scores against ESTIMATE and WES purity
- per-protein correlation with purity; compare R genes vs background
- stratify R genes by purity correlation and compare support
"""
import argparse
import csv
from pathlib import Path
import json
import os

import numpy as np
from scipy import stats

parser = argparse.ArgumentParser(description="Recompute the historical Figure S3 sensitivity.")
parser.add_argument("--inputs", type=Path, default=Path("inputs"))
parser.add_argument("--out", type=Path, default=Path("build/historical_purity"))
parser.add_argument("--membership", type=Path, default=Path(__file__).resolve().parents[1] / "results/audits/historical_purity_rna_membership.tsv")
args = parser.parse_args()
D, OUT = args.inputs, args.out
OUT.mkdir(parents=True, exist_ok=True)


def read_matrix(path):
    with open(path, encoding="utf-8", errors="replace") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        rows = [l.rstrip("\n").split("\t") for l in f]
    ids = np.array([r[0].split(".")[0] for r in rows])
    X = np.full((len(rows), len(hdr) - 1), np.nan)
    for i, r in enumerate(rows):
        for j, v in enumerate(r[1:]):
            if v and v != "NA":
                try:
                    X[i, j] = float(v)
                except ValueError:
                    pass
    return ids, np.array(hdr[1:]), X


pt, st, T = read_matrix(os.path.join(D, "cptac_pdac/PDAC_proteomics_gene_abundance_log2_reference_intensity_normalized_Tumor.txt"))
print("proteome:", T.shape)

ens2sym = {}
with open(os.path.join(D, "hgnc_complete_set.txt"), encoding="utf-8", errors="replace") as f:
    h = f.readline().rstrip("\n").split("\t")
    ie, isy = h.index("ensembl_gene_id"), h.index("symbol")
    for line in f:
        pp = line.rstrip("\n").split("\t")
        if len(pp) > max(ie, isy) and pp[ie]:
            ens2sym[pp[ie]] = pp[isy].upper()

sym = np.array([ens2sym.get(g, "") for g in pt])
keep = sym != ""
sym, T = sym[keep], T[keep, :]
# collapse duplicates by max mean
mean = np.nanmean(T, axis=1)
best = {}
for i, s in enumerate(sym):
    if s not in best or mean[i] > mean[best[s]]:
        best[s] = i
genes = np.array(list(best))
T = T[[best[g] for g in genes], :]
print("genes:", len(genes), "samples:", T.shape[1])

# phenotype
phen = {}
with open(os.path.join(D, "cptac_pdac/PDAC_phenotype.txt"), encoding="utf-8", errors="replace") as f:
    rd = csv.DictReader(f, delimiter="\t")
    cols = rd.fieldnames
    for r in rd:
        phen[r["idx"]] = r

samples = [s for s in st if s in phen]
cols_idx = [list(st).index(s) for s in samples]
Ts = T[:, cols_idx]
print("samples with phenotype:", len(samples))

purity = np.array([float(phen[s]["WES_purity"]) if phen[s]["WES_purity"] not in ("", "NA") else np.nan for s in samples])
stroma_est = np.array([float(phen[s]["ESTIMATE_StromalScore"]) if phen[s]["ESTIMATE_StromalScore"] not in ("", "NA") else np.nan for s in samples])

EPI = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "MUC1"]
STR = ["COL1A1", "COL1A2", "COL3A1", "COL5A1", "DCN", "LUM", "FAP", "ACTA2", "TAGLN", "PDGFRB", "THY1"]
gi = {g: i for i, g in enumerate(genes)}


def score(markers):
    idx = [gi[g] for g in markers if g in gi]
    Z = Ts[idx, :]
    Z = (Z - np.nanmean(Z, axis=1, keepdims=True)) / np.nanstd(Z, axis=1, keepdims=True)
    return np.nanmean(Z, axis=0), len(idx)


epi, ne = score(EPI)
strm, ns = score(STR)
print("markers used: epi %d, stroma %d" % (ne, ns))

mask = np.isfinite(purity) & np.isfinite(epi) & np.isfinite(strm)
print("with purity:", int(mask.sum()))
print("epi vs purity  rho=%.3f p=%.2g" % stats.spearmanr(epi[mask], purity[mask]))
print("stroma vs ESTIMATE rho=%.3f p=%.2g" % stats.spearmanr(strm[mask], stroma_est[mask]))
print("stroma vs purity rho=%.3f" % stats.spearmanr(strm[mask], purity[mask])[0])

# per-protein correlation with purity
rhos = np.full(len(genes), np.nan)
for i in range(len(genes)):
    x = Ts[i, mask]
    if np.isfinite(x).sum() > 30:
        rhos[i] = stats.spearmanr(x, purity[mask])[0]

with args.membership.open(encoding="utf-8-sig", newline="") as membership_stream:
    rows = list(csv.DictReader(membership_stream, delimiter="\t"))
R = {r["gene"] for r in rows if r["in_R"] == "1"}
R5 = {r["gene"] for r in rows if r["in_R_strict"] == "1"}
inset = np.array([g in R for g in genes])
inset5 = np.array([g in R5 for g in genes])

print("\nmedian |rho with purity|:")
print("  R      %.3f (n=%d)" % (np.nanmedian(np.abs(rhos[inset])), inset.sum()))
print("  R5     %.3f (n=%d)" % (np.nanmedian(np.abs(rhos[inset5])), inset5.sum()))
print("  others %.3f (n=%d)" % (np.nanmedian(np.abs(rhos[~inset & np.isfinite(rhos)])), (~inset).sum()))

# support by purity-correlation quartile among R
q = np.nanpercentile(np.abs(rhos[inset]), [25, 50, 75])
print("R |rho| quartiles:", np.round(q, 3))

res = dict(epi_vs_purity=float(stats.spearmanr(epi[mask], purity[mask])[0]),
           stroma_vs_estimate=float(stats.spearmanr(strm[mask], stroma_est[mask])[0]),
           stroma_vs_purity=float(stats.spearmanr(strm[mask], purity[mask])[0]),
           median_abs_rho_R=float(np.nanmedian(np.abs(rhos[inset]))),
           median_abs_rho_R5=float(np.nanmedian(np.abs(rhos[inset5]))),
           median_abs_rho_other=float(np.nanmedian(np.abs(rhos[~inset & np.isfinite(rhos)]))),
           n_purity=int(mask.sum()))
with open(os.path.join(OUT, "composition.json"), "w", encoding="utf-8") as f:
    json.dump(res, f, indent=1)

with (OUT / "TableS7_historical_purity.tsv").open("w", encoding="utf-8", newline="") as stream:
    writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
    writer.writerow(["group", "median_abs_rho_purity", "n"])
    writer.writerow(["R", f'{res["median_abs_rho_R"]:.3f}', int(inset.sum())])
    writer.writerow(["R_strict", f'{res["median_abs_rho_R5"]:.3f}', int(inset5.sum())])
    writer.writerow(["background", f'{res["median_abs_rho_other"]:.3f}', int((~inset).sum())])