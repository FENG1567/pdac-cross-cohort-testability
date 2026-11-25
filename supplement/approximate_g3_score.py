"""Model-level dedup demonstration on GSE62452.

Published coefficients (Chen et al., Table 5, doi:10.1093/g3journal/jkab296) are applied to GSE62452 expression;
Harrell C is computed on all tumours with survival, then on the subset that is NOT
one of the GSE28735 samples.
"""
import argparse
import csv
import gzip
from pathlib import Path
import json
import os
import re

import numpy as np

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--inputs", required=True, type=Path)
parser.add_argument("--map", required=True, type=Path)
parser.add_argument("--out", required=True, type=Path)
args = parser.parse_args()
IN, OUT = args.inputs, args.out
OUT.mkdir(parents=True, exist_ok=True)

COEF = {"FABP6": -0.0728, "REG4": 0.0570, "MAL": 0.0175, "KIF19": 0.0419,
        "CCK": -0.0377, "FER1L6": -0.0557, "SFTPA2": 0.0597, "GSTM1": 0.0331,
        "MUC17": -0.0076, "UPK1B": 0.0868, "HEPACAM2": -0.0877, "TNNT2": 0.0155,
        "NLRP2": -0.0920, "G6PC2": -0.0075, "KLK8": 0.0336}


def read_series(path):
    meta, cols, rows = {}, None, []
    started = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\n")
            if not started:
                if line.startswith("!series_matrix_table_begin"):
                    started = True
                    cols = [h.strip('"') for h in next(f).rstrip("\n").split("\t")[1:]]
                    continue
                if line.startswith("!"):
                    p = line.split("\t")
                    meta.setdefault(p[0], []).append([x.strip('"') for x in p[1:]])
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            rows.append(line.split("\t"))
    ids = [r[0].strip('"') for r in rows]
    X = np.array([[float(v) if v not in ("", "NA", "null") else np.nan for v in r[1:]] for r in rows])
    return meta, ids, cols, X


def load_annot(path):
    m, started = {}, False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not started:
                if line.startswith("ID\t"):
                    started = True
                continue
            if line.startswith("!platform_table_end"):
                break
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            sym = p[2].split("///")[0].strip().upper()
            if sym and sym != "---":
                m[p[0]] = sym
    return m


meta, ids, cols, X = read_series(os.path.join(IN, "GSE62452_series_matrix.txt.gz"))
annot = load_annot(os.path.join(IN, "GPL6244.annot.gz"))
rm = np.nanmean(X, axis=1)
best = {}
for i, pid in enumerate(ids):
    s = annot.get(pid)
    if s and (s not in best or rm[i] > rm[best[s]]):
        best[s] = i
genes = np.array([g for g in COEF if g in best])
print("model genes found:", genes.tolist(), "/", len(COEF))

# survival per sample
surv_months = meta["!Sample_characteristics_ch1"][3]
surv_status = meta["!Sample_characteristics_ch1"][4]
tissue = meta["!Sample_characteristics_ch1"][0]
titles = meta["!Sample_title"][0]
gsm = meta["!Sample_geo_accession"][0]

months, status = [], []
for a, b in zip(surv_months, surv_status):
    try:
        months.append(float(a.replace("survival months:", "").strip()))
    except ValueError:
        months.append(np.nan)
    v = b.replace("survival status:", "").strip()
    status.append(1.0 if v == "1" else (0.0 if v == "0" else np.nan))
months, status = np.array(months), np.array(status)

is_tumor = np.array(["tumor" in t.lower() and "nontumor" not in t.lower() for t in tissue]) | \
           np.array([bool(re.search(r"T\w*$", t)) and "Tn" not in t for t in titles])
print("tumours:", int(is_tumor.sum()), "with survival:", int((is_tumor & np.isfinite(months) & np.isfinite(status)).sum()))

expr = X[[best[g] for g in genes], :]
pi = np.array([COEF[g] for g in genes]) @ expr

# overlap: GSE62452 columns 1..90 correspond to GSE28735 samples
mapping = {}
with args.map.open(encoding="utf-8", newline="") as f:
    for r in csv.DictReader(f, delimiter="\t"):
        mapping[r["gse62452_gsm"]] = r

overlap_gsm = set(mapping)
is_overlap = np.array([g in overlap_gsm for g in gsm])

use = is_tumor & np.isfinite(months) & np.isfinite(status) & np.isfinite(pi)
print("usable tumours:", int(use.sum()), "of which overlap:", int((use & is_overlap).sum()))


def harrell_c(score, time, event):
    n = len(score)
    num = den = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            if time[i] == time[j]:
                continue
            if time[i] < time[j] and event[i] == 1:
                den += 1
                num += 1 if score[i] > score[j] else (0.5 if score[i] == score[j] else 0)
            elif time[j] < time[i] and event[j] == 1:
                den += 1
                num += 1 if score[j] > score[i] else (0.5 if score[j] == score[i] else 0)
    return num / den if den else np.nan


def boot_c(score, time, event, n=500, seed=1):
    rng = np.random.default_rng(seed)
    vals = []
    n_ = len(score)
    for _ in range(n):
        idx = rng.integers(0, n_, n_)
        vals.append(harrell_c(score[idx], time[idx], event[idx]))
    return np.nanpercentile(vals, [2.5, 97.5])


res = {}
for label, mask in (("all", use), ("overlap_only", use & is_overlap), ("unique_only", use & ~is_overlap)):
    s, t, e = pi[mask], months[mask], status[mask]
    c = harrell_c(s, t, e)
    ci = boot_c(s, t, e)
    res[label] = dict(n=int(mask.sum()), c=float(c), lo=float(ci[0]), hi=float(ci[1]))
    print(label, res[label])

with open(os.path.join(OUT, "dedup_cindex.json"), "w", encoding="utf-8") as f:
    json.dump(res, f, indent=1)
