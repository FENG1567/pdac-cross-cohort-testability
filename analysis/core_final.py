#!/usr/bin/env python3
"""Reproducible four-cohort RNA-to-protein support and ranking analysis.

This standalone runner reads inputs, preserves the frozen RNA sets,
deduplicates protein symbols by highest mean abundance, writes early RNA-only
top-500 ranking lists, then builds all point summaries and patient bootstraps.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
import re
import sys
import time
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


RULES = ("single_A", "single_B", "equal_weight", "continuous", "double_sig")
KS = (100, 250, 500, 1000, 2000)
MIN_PAIRS = 10
PRIMARY_Q = 0.05
SEED = 20261123
BOOTSTRAPS = 1000


def log(*parts):
    stamp = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    print(" ".join([stamp] + [str(x) for x in parts]), flush=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def parse_num(v):
    if v is None:
        return np.nan
    if isinstance(v, (int, float, np.number)):
        return float(v) if np.isfinite(v) else np.nan
    s = str(v).strip()
    if not s or s.upper() in {"NA", "NAN", "NULL", "NONE"}:
        return np.nan
    try:
        x = float(s)
        return x if np.isfinite(x) else np.nan
    except (TypeError, ValueError):
        return np.nan


def clean_symbol(s):
    if s is None:
        return ""
    s = str(s).strip().upper()
    if s in {"", "NA", "NAN", "NULL", "---", "-"}:
        return ""
    return s


def read_series_matrix(path: Path):
    meta = {}
    cols, rows = None, []
    started = False
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.rstrip("\r\n")
            if not started:
                if line.startswith("!series_matrix_table_begin"):
                    started = True
                    cols = [x.strip('"') for x in next(f).rstrip("\r\n").split("\t")[1:]]
                    continue
                if line.startswith("!"):
                    p = line.split("\t")
                    meta.setdefault(p[0], []).append([x.strip('"') for x in p[1:]])
                continue
            if line.startswith("!series_matrix_table_end"):
                break
            rows.append(line.split("\t"))
    if cols is None:
        raise ValueError(f"No series_matrix_table_begin in {path}")
    ids = [r[0].strip('"') for r in rows]
    X = np.full((len(rows), len(cols)), np.nan, float)
    for i, row in enumerate(rows):
        for j, v in enumerate(row[1:len(cols) + 1]):
            X[i, j] = parse_num(v.strip('"'))
    return meta, ids, cols, X


def load_annotation(path: Path):
    out = {}
    started = False
    with gzip.open(path, "rt", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not started:
                if line.startswith("ID\t"):
                    started = True
                continue
            if line.startswith("!platform_table_end"):
                break
            p = line.rstrip("\r\n").split("\t")
            if len(p) < 3:
                continue
            syms = [clean_symbol(s) for s in p[2].strip().split("///")]
            syms = [s for s in syms if s]
            if syms:
                out[p[0]] = syms
    return out


def collapse_probes(X, ids, annot):
    rowmean = np.divide(np.nansum(X, axis=1), np.sum(np.isfinite(X), axis=1),
                        out=np.full(X.shape[0], np.nan), where=np.sum(np.isfinite(X), axis=1) > 0)
    best = {}
    for i, pid in enumerate(ids):
        for gene in annot.get(pid, []):
            old = best.get(gene)
            if old is None or (np.isfinite(rowmean[i]) and
                               (not np.isfinite(rowmean[old]) or rowmean[i] > rowmean[old])):
                best[gene] = i
    genes = np.array(sorted(best), dtype=object)
    take = np.array([best[g] for g in genes], dtype=int)
    return genes, X[take, :]


def make_rna_cohort(tag, series_path, annot_path, exclude_gse15471_qc=False):
    meta, ids, cols, X = read_series_matrix(series_path)
    titles = meta["!Sample_title"][0]
    if tag == "A_GSE28735":
        pat = {}
        for i, title in enumerate(titles):
            m = re.search(r"patient sample (\d+)", title, flags=re.I)
            if not m:
                continue
            pat.setdefault(int(m.group(1)), {})["N" if "nontumor" in title.lower() else "T"] = i
        pairs = [(k, v["T"], v["N"]) for k, v in pat.items() if "T" in v and "N" in v]
        pairs.sort(key=lambda x: x[0])
    else:
        pat = {}
        for i, title in enumerate(titles):
            m = re.match(r"^([TN])(\d+)", title)
            if not m or "_rep" in title.lower():
                continue
            pat.setdefault(m.group(2), {})[m.group(1)] = i
        pairs = [(k, v["T"], v["N"]) for k, v in pat.items() if "T" in v and "N" in v]
        pairs.sort(key=lambda x: int(x[0]))
        rawmax = float(np.nanmax(X))
        if rawmax > 50:
            X = np.log2(np.clip(X, 0, None) + 1.0)
    if exclude_gse15471_qc:
        assert tag == "B_GSE15471" and len(pairs) == 36
        diagnostic = {"GSM388111", "GSM388150"}
        retained = [pair for pair in pairs if not ({cols[pair[1]], cols[pair[2]]} & diagnostic)]
        removed = [pair for pair in pairs if {cols[pair[1]], cols[pair[2]]} == diagnostic]
        if len(removed) != 1 or len(retained) != 35:
            raise ValueError("Reported GSE15471 quality-control samples do not form a single retained pair")
        # Keep excluded arrays out of the expression-based probe selection as well.
        keep = np.array([index for _, tumour, normal in retained for index in (tumour, normal)])
        X = X[:, keep]
        pairs = [(pair[0], 2*i, 2*i + 1) for i, pair in enumerate(retained)]
        log("GSE15471 quality-control pair confirmed", f"retained_pairs={len(pairs)}")
    genes, Xg = collapse_probes(X, ids, load_annotation(annot_path))
    ti = np.array([x[1] for x in pairs], dtype=int)
    ni = np.array([x[2] for x in pairs], dtype=int)
    D = Xg[:, ti] - Xg[:, ni]
    return dict(name=tag, genes=genes, X=Xg, Dm=D, pair_names=[str(x[0]) for x in pairs],
                pairs=len(pairs), ti=ti, ni=ni)


def bh(p):
    p = np.asarray(p, float)
    q = np.full(p.shape, np.nan, float)
    ok = np.isfinite(p)
    pv = p[ok]
    n = len(pv)
    if n == 0:
        return q
    order = np.argsort(pv, kind="mergesort")
    vals = pv[order] * n / np.arange(1, n + 1, dtype=float)
    vals = np.minimum.accumulate(vals[::-1])[::-1]
    vals = np.clip(vals, 0, 1)
    tmp = np.empty(n, float)
    tmp[order] = vals
    q[ok] = tmp
    return q


def paired_stats(Dm, min_pairs=MIN_PAIRS):
    Dm = np.asarray(Dm, float)
    if Dm.ndim != 2:
        raise ValueError("delta matrix must be two-dimensional")
    valid = np.isfinite(Dm)
    n = valid.sum(axis=1).astype(int)
    sm = np.where(valid, Dm, 0.0).sum(axis=1)
    mean = np.divide(sm, n, out=np.full(Dm.shape[0], np.nan), where=n > 0)
    centered = np.where(valid, Dm - mean[:, None], 0.0)
    ss = np.square(centered).sum(axis=1)
    var = np.divide(ss, n - 1, out=np.full(Dm.shape[0], np.nan), where=n > 1)
    sd = np.sqrt(np.maximum(var, 0.0))
    se = np.divide(sd, np.sqrt(n), out=np.full(Dm.shape[0], np.nan), where=n > 0)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        t = mean / se
        p = 2.0 * stats.t.sf(np.abs(t), np.maximum(n - 1, 1))
    zero_sd = (sd == 0) & (n >= 2) & np.isfinite(mean)
    p[zero_sd & (mean != 0)] = 0.0
    p[zero_sd & (mean == 0)] = np.nan
    bad = (n < min_pairs) | ~np.isfinite(p)
    raw_mean = mean.copy()
    p[bad] = np.nan
    mean[bad] = np.nan
    se[bad] = np.nan
    q = bh(p)
    return dict(mean=mean, raw_mean=raw_mean, p=p, q=q, n=n, sd=sd, se=se)


def scores_and_orders(genes, dA, seA, qA, dB, seB, qB):
    with np.errstate(divide="ignore", invalid="ignore"):
        wA = 1.0 / np.square(np.clip(seA, 1e-6, None))
        wB = 1.0 / np.square(np.clip(seB, 1e-6, None))
        comp = (wA * dA + wB * dB) / (wA + wB)
    in_double = (qA < PRIMARY_Q) & (qB < PRIMARY_Q) & (np.sign(dA) == np.sign(dB))
    data = {
        "single_A": (np.abs(dA), np.sign(dA)),
        "single_B": (np.abs(dB), np.sign(dB)),
        "equal_weight": (np.abs(dA + dB), np.sign(dA + dB)),
        "continuous": (np.abs(comp), np.sign(comp)),
        "double_sig": (np.where(in_double, 1e6, 0.0) + np.abs(comp), np.sign(comp)),
    }
    orders = {}
    for rule in RULES:
        score, direction = data[rule]
        clean = np.where(np.isfinite(score), score, -np.inf)
        # Lexicographic secondary key is gene symbol ascending; deterministic at ties.
        order = np.lexsort((np.asarray(genes, dtype=str), -clean))
        orders[rule] = dict(score=np.asarray(score), direction=np.asarray(direction), order=order)
    return orders


def read_frozen(path):
    df = pd.read_csv(path, sep="\t", dtype={"gene": str})
    if not {"gene", "dir"}.issubset(df.columns):
        raise ValueError(f"frozen set lacks gene/dir columns: {path}")
    genes = df["gene"].astype(str).map(clean_symbol)
    if genes.duplicated().any():
        raise ValueError(f"duplicate gene symbols in frozen set {path}")
    return {g: int(d) for g, d in zip(genes, df["dir"].astype(int)) if g}


def write_tsv(df, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, na_rep="NA", float_format="%.10g",
              encoding="utf-8", lineterminator="\n")


def mean_abundance(X):
    n = np.isfinite(X).sum(axis=1)
    s = np.where(np.isfinite(X), X, 0.0).sum(axis=1)
    return np.divide(s, n, out=np.full(X.shape[0], np.nan), where=n > 0)


def select_highest_mean(cohort, raw_symbols, abundance_X, delta_X, source_labels):
    """Select highest mean abundance row; source order breaks exact ties (first wins)."""
    symbols = [clean_symbol(s) for s in raw_symbols]
    if len(symbols) != abundance_X.shape[0] or delta_X.shape[0] != len(symbols):
        raise ValueError(f"row count mismatch in {cohort}")
    avg = mean_abundance(np.asarray(abundance_X, float))
    groups = defaultdict(list)
    blank_rows = 0
    for i, g in enumerate(symbols):
        if g:
            groups[g].append(i)
        else:
            blank_rows += 1
    selected, last = {}, {}
    details = []
    tied_max = 0
    all_nan = 0
    for gene, ix in groups.items():
        vals = avg[ix]
        if not np.isfinite(vals).any():
            best_i = ix[0]
            all_nan += 1
        else:
            vmax = np.nanmax(vals)
            tied = [i for i in ix if np.isfinite(avg[i]) and np.isclose(avg[i], vmax, rtol=0, atol=0)]
            if len(tied) > 1:
                tied_max += 1
            best_i = min(tied)
        last_i = ix[-1]
        selected[gene] = best_i
        last[gene] = last_i
        if len(ix) > 1:
            details.append(dict(
                cohort=cohort, gene=gene, source_rows=len(ix),
                selected_row_index0=best_i, selected_row_label=source_labels[best_i],
                selected_mean_abundance=avg[best_i], last_row_index0=last_i,
                last_row_label=source_labels[last_i], last_mean_abundance=avg[last_i],
                selected_is_last_row=bool(best_i == last_i),
                selection_rule="max finite mean abundance; exact ties -> earliest source row"))
    genes = np.array(sorted(selected), dtype=object)
    take = np.array([selected[g] for g in genes], dtype=int)
    take_last = np.array([last[g] for g in genes], dtype=int)
    summary = dict(cohort=cohort, rows_with_symbol=len(symbols) - blank_rows,
                   blank_or_invalid_symbol_rows=blank_rows, unique_symbols_before=len(groups),
                   duplicate_symbols=sum(len(v) > 1 for v in groups.values()),
                   rows_removed_by_dedup=sum(max(0, len(v) - 1) for v in groups.values()),
                   rows_before_dedup=len(symbols) - blank_rows, rows_after_dedup=len(genes),
                   unique_symbols_after=len(genes), all_nan_groups=all_nan,
                   groups_tied_at_max_mean=tied_max,
                   historical_lastrow_groups_whose_selected_row_differs=int(sum(
                       selected[g] != last[g] for g in selected)),
                   selection_rule="max mean abundance across observed samples; exact ties -> earliest source row")
    return dict(genes=genes, Dm=np.asarray(delta_X, float)[take, :],
                abundance=avg[take], source_labels=[source_labels[i] for i in take],
                selected_source_row_index=take, last_Dm=np.asarray(delta_X, float)[take_last, :],
                last_abundance=avg[take_last], summary=summary, details=details)


def read_cct(path):
    with path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.reader(f, delimiter="\t")
        header = next(reader)
        rows = list(reader)
    if not header:
        raise ValueError(f"empty CCT {path}")
    samples = [str(x).strip() for x in header[1:]]
    ids, vals = [], np.full((len(rows), len(samples)), np.nan, float)
    for i, row in enumerate(rows):
        if not row:
            ids.append("")
            continue
        ids.append(str(row[0]).split(".")[0].strip())
        for j, v in enumerate(row[1:len(samples) + 1]):
            vals[i, j] = parse_num(v)
    return ids, samples, vals


def load_fudan(path):
    from openpyxl import load_workbook
    wb = load_workbook(path, read_only=True, data_only=True)
    if "Table S3A " in wb.sheetnames:
        ws = wb["Table S3A "]
    elif "Table S3A" in wb.sheetnames:
        ws = wb["Table S3A"]
    else:
        raise ValueError(f"Fudan sheet missing; observed {wb.sheetnames}")
    it = ws.iter_rows(values_only=True)
    hdr = [str(x).strip() if x is not None else "" for x in next(it)]
    keep = [i for i, x in enumerate(hdr) if i > 0 and x]
    samples = [hdr[i] for i in keep]
    symbols, rows = [], []
    for row in it:
        symbols.append(clean_symbol(row[0] if row else ""))
        rows.append([parse_num(row[i] if i < len(row) else None) for i in keep])
    wb.close()
    raw = np.asarray(rows, float)
    X = np.log2(raw + 1.0)
    tcol = {c[:-1]: i for i, c in enumerate(samples) if c.endswith("T")}
    ncol = {c[:-1]: i for i, c in enumerate(samples) if c.endswith("N")}
    pairs = sorted(set(tcol) & set(ncol))
    if not pairs:
        raise ValueError("Fudan paired T/N columns were not found")
    ti, ni = [tcol[p] for p in pairs], [ncol[p] for p in pairs]
    return dict(name="Fudan", raw_symbols=symbols, abundance_X=X, Dm=X[:, ti] - X[:, ni],
                pair_names=pairs, source_labels=[f"xlsx_data_row_{i + 2}" for i in range(len(symbols))])


def load_cptac(tumor_path, normal_path):
    tid, ts, T = read_cct(tumor_path)
    nid, ns, N = read_cct(normal_path)
    pairs = sorted(set(ts) & set(ns))
    if not pairs:
        raise ValueError("CPTAC paired tumor/normal sample ids were not found")
    ti = [ts.index(p) for p in pairs]
    ni = [ns.index(p) for p in pairs]
    tg, ng = defaultdict(list), defaultdict(list)
    for i, g in enumerate(tid):
        g = clean_symbol(g)
        if g:
            tg[g].append(i)
    for i, g in enumerate(nid):
        g = clean_symbol(g)
        if g:
            ng[g].append(i)
    syms, abund_rows, delta_rows, labels = [], [], [], []
    for gene in sorted(set(tg) & set(ng)):
        n = max(len(tg[gene]), len(ng[gene]))
        for j in range(n):
            tr = tg[gene][j] if j < len(tg[gene]) else None
            nr = ng[gene][j] if j < len(ng[gene]) else None
            ab = np.concatenate([T[tr, :] if tr is not None else np.full(T.shape[1], np.nan),
                                 N[nr, :] if nr is not None else np.full(N.shape[1], np.nan)])
            td = T[tr, ti] if tr is not None else np.full(len(pairs), np.nan)
            nd = N[nr, ni] if nr is not None else np.full(len(pairs), np.nan)
            syms.append(gene)
            abund_rows.append(ab)
            delta_rows.append(td - nd)
            labels.append(f"tumor_row_{tr + 2 if tr is not None else 'NA'}|normal_row_{nr + 2 if nr is not None else 'NA'}")
    return dict(name="CPTAC", raw_symbols=syms, abundance_X=np.asarray(abund_rows, float),
                Dm=np.asarray(delta_rows, float), pair_names=pairs, source_labels=labels)


def load_apgi(matrix_path, hgnc_path):
    with matrix_path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.reader(f)
        hdr = next(reader)
        cut = hdr.index("Prep_LIMS_ID")
        histology_col = hdr.index("Histology")
        uids = [x.strip() for x in hdr[1:cut]]
        samples, rows, sample_histology = [], [], []
        for row in reader:
            if not row:
                continue
            samples.append(row[0].strip())
            rows.append([parse_num(v) for v in row[1:cut]])
            sample_histology.append(row[histology_col].strip() if len(row) > histology_col else "")
    X = np.asarray(rows, float).T
    if len(uids) != X.shape[0]:
        raise ValueError("APGI UID/sample matrix dimensions disagree")
    uni2sym = {}
    with hgnc_path.open("r", encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            sym = clean_symbol(row.get("symbol", ""))
            if not sym:
                continue
            for uid in str(row.get("uniprot_ids", "") or "").split("|"):
                uid = uid.strip()
                if uid:
                    uni2sym.setdefault(uid, sym)
    pat = {}
    for i, s in enumerate(samples):
        m = re.match(r"^(\d+)_(Tumour|Normal)$", s)
        if m:
            pat.setdefault(m.group(1), {})["T" if m.group(2) == "Tumour" else "N"] = i
    patients = sorted((k for k, v in pat.items() if "T" in v and "N" in v), key=int)
    if not patients:
        raise ValueError("APGI paired Tumour/Normal samples were not found")
    ti = [pat[p]["T"] for p in patients]
    ni = [pat[p]["N"] for p in patients]
    pair_histology = {p: sample_histology[pat[p]["T"]] for p in patients}
    strict_patients = [p for p in patients if pair_histology[p] == "Pancreatic Ductal Adenocarcinoma"]
    eligibility = [dict(patient_id=p, tumor_sample=samples[pat[p]["T"]], normal_sample=samples[pat[p]["N"]],
                        tumor_histology=pair_histology[p], strict_pdac=pair_histology[p] == "Pancreatic Ductal Adenocarcinoma",
                        inclusion_rule="tumor-row Histology exactly equals Pancreatic Ductal Adenocarcinoma")
                   for p in patients]
    symbols = [uni2sym.get(uid.split("-")[0], "") for uid in uids]
    return dict(name="APGI", raw_symbols=symbols, abundance_X=X, Dm=X[:, ti] - X[:, ni],
                pair_names=patients, source_labels=uids, pair_histology=pair_histology,
                strict_patient_indices=np.array([patients.index(p) for p in strict_patients], dtype=int),
                eligibility=eligibility)


def load_ruijin(path):
    # This workbook has a damaged worksheet dimension/drawing reference. In
    # read_only mode openpyxl can report a 1x1 sheet, so use the non-streaming
    # pandas/openpyxl path validated on this input (4,787 x 282 including genes).
    df = pd.read_excel(path, engine="openpyxl")
    cols = [str(x).strip() if x is not None else "" for x in df.columns]
    symbols = [clean_symbol(x) for x in df.iloc[:, 0].tolist()]
    X = df.iloc[:, 1:].apply(pd.to_numeric, errors="coerce").to_numpy(dtype=float)
    labels = [f"xlsx_row_{i + 2}" for i in range(len(symbols))]
    colnames = cols[1:]
    pat = {}
    for i, col in enumerate(colnames):
        m = re.match(r"^(RJ-\d+-\d+)_([TN])_", col)
        if m:
            pat.setdefault(m.group(1), {})[m.group(2)] = i
    patients = sorted((k for k, v in pat.items() if "T" in v and "N" in v))
    if not patients:
        raise ValueError("Ruijin paired T/N sample columns were not found")
    ti = [pat[p]["T"] for p in patients]
    ni = [pat[p]["N"] for p in patients]
    return dict(name="Ruijin", raw_symbols=symbols, abundance_X=X, Dm=X[:, ti] - X[:, ni],
                pair_names=patients, source_labels=labels)


def bootstrap_stats(cohort, indices):
    return paired_stats(cohort["Dm"][:, indices], min_pairs=MIN_PAIRS)


def align_stats_to_master(genes, master, stat):
    pos = {g: i for i, g in enumerate(genes)}
    take = np.array([pos.get(g, -1) for g in master], int)
    out = {}
    for k in ("mean", "p", "q", "n"):
        dtype = float if k != "n" else float
        v = np.full(len(master), np.nan, dtype=dtype)
        ok = take >= 0
        v[ok] = stat[k][take[ok]]
        out[k] = v
    out["has_protein"] = take >= 0
    return out


def make_manifest(inputs, output, only_files=None):
    entries = []
    paths = ([inputs / rel for rel in only_files] if only_files is not None else sorted(inputs.rglob("*")))
    for p in sorted(paths):
        if p.is_file():
            entries.append(dict(file=str(p.relative_to(inputs)), bytes=p.stat().st_size, sha256=sha256(p)))
    write_tsv(pd.DataFrame(entries), output / "input_manifest.tsv")
    return entries


def export_early_rankings(master, dA, seA, qA, dB, seB, qB, frozen_R, outdir, input_hashes):
    orders = scores_and_orders(master, dA, seA, qA, dB, seB, qB)
    records = []
    for rule in RULES:
        info = orders[rule]
        for rank, idx in enumerate(info["order"][:500], start=1):
            records.append(dict(rule=rule, rank=rank, gene=str(master[idx]), score=info["score"][idx],
                                direction=int(np.sign(info["direction"][idx])),
                                in_frozen_R=bool(master[idx] in frozen_R)))
    df = pd.DataFrame(records)
    outdir.mkdir(parents=True, exist_ok=True)
    canonical = outdir / "unrestricted_top500_all_methods.tsv"
    write_tsv(df, canonical)
    bio = outdir / "biotext_top500_all_methods.tsv"
    prm = outdir / "prm_top500_all_methods.tsv"
    df.to_csv(bio, sep="\t", index=False, na_rep="NA", float_format="%.10g", encoding="utf-8", lineterminator="\n")
    df.to_csv(prm, sep="\t", index=False, na_rep="NA", float_format="%.10g", encoding="utf-8", lineterminator="\n")
    all_scores = []
    all_ranks = {}
    for rule in RULES:
        rank = np.empty(len(master), dtype=int)
        rank[orders[rule]["order"]] = np.arange(1, len(master) + 1)
        all_ranks[rule] = rank
    for i, gene in enumerate(master):
        row = dict(gene=str(gene), d_A=dA[i], se_A=seA[i], q_A=qA[i], d_B=dB[i], se_B=seB[i], q_B=qB[i],
                   in_frozen_R=gene in frozen_R)
        for rule in RULES:
            info = orders[rule]
            row[f"score_{rule}"] = info["score"][i]
            row[f"direction_{rule}"] = int(np.sign(info["direction"][i]))
            # rank is added below from the deterministic order.
        for rule in RULES:
            row[f"rank_{rule}"] = int(all_ranks[rule][i])
        all_scores.append(row)
    write_tsv(pd.DataFrame(all_scores), outdir / "rna_scores_and_all_ranks.tsv")
    meta = dict(master_genes=int(len(master)), frozen_R_genes=int(len(frozen_R)),
                rules=list(RULES), tie_break="score descending, gene symbol ascending",
                rank_universe="all shared A/B RNA genes; unrestricted by protein measurement or significance",
                cutoff=500, formulas={
                    "single_A": "abs(d_A)", "single_B": "abs(d_B)",
                    "equal_weight": "abs(d_A+d_B)",
                    "continuous": "abs(inverse-variance weighted mean(d_A,d_B)); SE from paired deltas",
                    "double_sig": "1e6 priority for frozen protocol-R criterion plus abs(continuous score)"},
                sha256={p.name: sha256(p) for p in (canonical, bio, prm)},
                frozen_protocol_sha256=input_hashes.get("R_protocol.tsv"))
    (outdir / "ranking_manifest.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    return orders, df


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", required=True, type=Path)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--rankings-only", action="store_true",
                    help="write frozen unrestricted RNA ranking artifacts, then exit before protein reads")
    ap.add_argument("--cptac-histology", type=Path,
                    help="Complete case-key audit (case_key, histology, keep_pdac) for paired CPTAC cases")
    ap.add_argument("--historical-sets", type=Path, help="Historical R sets used for membership comparison")
    args = ap.parse_args()
    inputs = args.inputs.resolve()
    outdir = args.out.resolve()
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "logs").mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    log("START four-cohort core run", f"inputs={inputs}", f"out={outdir}",
        f"bootstrap={BOOTSTRAPS}", f"seed={SEED}", f"min_pairs={MIN_PAIRS}")

    ranking_inputs = [Path("GSE28735_series_matrix.txt.gz"), Path("GSE15471_series_matrix.txt.gz"),
                      Path("GPL6244.annot.gz"), Path("GPL570.annot.gz")]
    protein_inputs = [Path("13045_2022_1384_MOESM25_ESM.xlsx"),
                      Path("cptac_orig/proteomics_gene_level_MD_abundance_tumor.cct"),
                      Path("cptac_orig/proteomics_gene_level_MD_abundance_normal.cct"),
                      Path("apgi/protein_matrix.csv"), Path("hgnc_complete_set.txt"),
                      Path("node_ruijin/20221227_PDAC_PRO_exp.xlsx")]
    selected_inputs = ranking_inputs if args.rankings_only else ranking_inputs + protein_inputs
    manifest = make_manifest(inputs, outdir, selected_inputs)
    input_hashes = {Path(x["file"]).name: x["sha256"] for x in manifest}
    historical_protocol = args.historical_sets / "R_protocol.tsv" if args.historical_sets else None
    historical_strict = args.historical_sets / "R_strict.tsv" if args.historical_sets else None
    historical_comparison_available = bool(historical_protocol and historical_strict and
                                       historical_protocol.is_file() and historical_strict.is_file())
    if args.historical_sets and not historical_comparison_available:
        raise FileNotFoundError("--historical-sets must contain R_protocol.tsv and R_strict.tsv")
    if historical_comparison_available:
        R = read_frozen(historical_protocol)
        R5 = read_frozen(historical_strict)
        if not set(R5).issubset(R):
            raise ValueError("frozen R_strict is not a subset of R_protocol")
        log("optional historical RNA comparators", f"R={len(R)}", f"R_strict={len(R5)}",
            f"R_sha256={sha256(historical_protocol)}", f"R_strict_sha256={sha256(historical_strict)}")
    else:
        # The primary R/R_strict sets are recalculated below from the paired
        # RNA cohorts. Historical sets are used only for optional membership comparison.
        R, R5 = {}, {}
        log("no historical RNA comparator supplied; recalculating primary sets from RNA cohorts")

    A = make_rna_cohort("A_GSE28735", inputs / "GSE28735_series_matrix.txt.gz", inputs / "GPL6244.annot.gz")
    B = make_rna_cohort("B_GSE15471", inputs / "GSE15471_series_matrix.txt.gz", inputs / "GPL570.annot.gz",
                        exclude_gse15471_qc=True)
    statA = paired_stats(A["Dm"])
    statB = paired_stats(B["Dm"])
    rna_qc_rows = []
    for tag, cohort, stat in (("A_GSE28735", A, statA), ("B_GSE15471", B, statB)):
        rna_qc_rows.append(dict(cohort=tag, genes=len(cohort["genes"]), paired_patients=cohort["pairs"],
            genes_with_any_missing_delta=int(np.sum(np.any(~np.isfinite(cohort["Dm"]), axis=1))),
            complete_pair_min=MIN_PAIRS, estimable_BH_tests=int(np.isfinite(stat["p"]).sum()),
            zero_sd_nonzero=int(np.sum((stat["n"] >= MIN_PAIRS) & (stat["sd"] == 0) &
                                       np.isfinite(stat["mean"]) & (stat["mean"] != 0))),
            zero_sd_zero_mean=int(np.sum((stat["n"] >= MIN_PAIRS) & (stat["sd"] == 0) &
                                         np.isfinite(stat["raw_mean"]) & (stat["raw_mean"] == 0))),
            testable_rule="paired log-scale/no-imputation; finite p after min-pair rule enters cohort-wide BH"))
    write_tsv(pd.DataFrame(rna_qc_rows), outdir / "rna_bh_family_summary.tsv")
    ia, ib = {g: i for i, g in enumerate(A["genes"])}, {g: i for i, g in enumerate(B["genes"])}
    master = np.array(sorted(set(ia) & set(ib)), dtype=object)
    if len(master) != 19187:
        log("WARNING expected historical shared RNA count=19187, observed", len(master))
    selA = np.array([ia[g] for g in master], int)
    selB = np.array([ib[g] for g in master], int)
    dA, seA, qA = statA["mean"][selA], statA["se"][selA], statA["q"][selA]
    dB, seB, qB = statB["mean"][selB], statB["se"][selB], statB["q"][selB]
    calcR = (qA < PRIMARY_Q) & (qB < PRIMARY_Q) & (np.sign(dA) == np.sign(dB))
    calcRdict = {str(g): int(np.sign(dA[i] + dB[i])) for i, g in enumerate(master) if calcR[i]}
    calcR5 = calcR & (np.abs(dA) >= 0.5) & (np.abs(dB) >= 0.5)
    calcR5dict = {str(g): int(np.sign(dA[i] + dB[i])) for i, g in enumerate(master) if calcR5[i]}
    if historical_comparison_available:
        frozen_compare = dict(historical_comparison_available=True,
                              frozen_R=len(R), calculated_R=int(calcR.sum()),
                              frozen_not_calculated=len(set(R) - set(calcRdict)),
                              calculated_not_frozen=len(set(calcRdict) - set(R)),
                              direction_mismatches=sum(R.get(g) != calcRdict.get(g) for g in set(R) & set(calcRdict)),
                              frozen_R_strict=len(R5), calculated_R_strict=int(calcR5.sum()),
                              frozen_strict_not_calculated=len(set(R5) - set(calcR5dict)),
                              calculated_strict_not_frozen=len(set(calcR5dict) - set(R5)),
                              strict_direction_mismatches=sum(R5.get(g) != calcR5dict.get(g) for g in set(R5) & set(calcR5dict)))
    else:
        frozen_compare = dict(historical_comparison_available=False, frozen_R=None,
                              calculated_R=int(calcR.sum()), frozen_not_calculated=None,
                              calculated_not_frozen=None, direction_mismatches=None,
                              frozen_R_strict=None, calculated_R_strict=int(calcR5.sum()),
                              frozen_strict_not_calculated=None, calculated_strict_not_frozen=None,
                              strict_direction_mismatches=None)
    R, R5 = calcRdict, calcR5dict
    frozen_compare["primary_definition"] = "GSE15471 documented QC pair omitted; R/R_strict recalculated from paired RNA cohorts"
    frozen_compare["historical_comparison"] = ("performed" if historical_comparison_available else "not supplied")
    (outdir / "frozen_set_validation.json").write_text(json.dumps(frozen_compare, indent=2), encoding="utf-8")
    log("RNA loaded", f"A_pairs={A['pairs']}", f"A_genes={len(A['genes'])}",
        f"B_pairs={B['pairs']}", f"B_genes={len(B['genes'])}", f"common={len(master)}",
        f"frozen_R_match={frozen_compare}")

    # Export unrestricted rankings before any protein input is read or evaluated.
    ranking_dir = outdir.parent / "rankings"
    orders, rank500_df = export_early_rankings(master, dA, seA, qA, dB, seB, qB, set(R), ranking_dir, input_hashes)
    log("[RANKINGS_READY]", ranking_dir / "biotext_top500_all_methods.tsv",
        ranking_dir / "prm_top500_all_methods.tsv")
    rank_frame = pd.read_csv(ranking_dir / "rna_scores_and_all_ranks.tsv", sep="\t", dtype={"gene": str})
    rank_frame["R_direction"] = rank_frame["gene"].map(R).fillna(0).astype(int)
    rank_frame["R_strict_direction"] = rank_frame["gene"].map(R5).fillna(0).astype(int)
    write_tsv(rank_frame, outdir / "rna_scores_all_ranks.tsv")
    write_tsv(pd.DataFrame([dict(gene=g, direction=int(d)) for g, d in sorted(R.items())]), outdir / "R_protocol.tsv")
    write_tsv(pd.DataFrame([dict(gene=g, direction=int(d)) for g, d in sorted(R5.items())]), outdir / "R_strict.tsv")
    if args.rankings_only:
        early = dict(completed_utc=datetime.now(timezone.utc).isoformat(), status="rankings_only_complete",
                     seed=SEED, shared_RNA_gene_count=len(master), frozen_protocol_R_count=len(R),
                     ranking_files={p.name: sha256(p) for p in (ranking_dir / "biotext_top500_all_methods.tsv",
                         ranking_dir / "prm_top500_all_methods.tsv", ranking_dir / "rna_scores_and_all_ranks.tsv")})
        (outdir / "early_ranking_run.json").write_text(json.dumps(early, indent=2), encoding="utf-8")
        log("DONE rankings-only", f"elapsed_s={round(time.time() - t0, 1)}")
        return

    raw = []
    raw.append(load_fudan(inputs / "13045_2022_1384_MOESM25_ESM.xlsx"))
    cptac = load_cptac(inputs / "cptac_orig" / "proteomics_gene_level_MD_abundance_tumor.cct",
                        inputs / "cptac_orig" / "proteomics_gene_level_MD_abundance_normal.cct")
    if not args.cptac_histology:
        raise ValueError("Full analysis requires the verified CPTAC histology key table")
    clinical = pd.read_csv(args.cptac_histology, sep="\t", dtype={"case_key": str})
    assert clinical.case_key.is_unique and set(clinical.case_key) == set(cptac["pair_names"])
    clinical_keys = set(clinical.loc[clinical.keep_pdac == 1, "case_key"])
    assert len(clinical_keys) >= 60
    kept = [i for i, case in enumerate(cptac["pair_names"]) if case in clinical_keys]
    cptac["Dm"] = cptac["Dm"][:, kept]
    cptac["pair_names"] = [cptac["pair_names"][i] for i in kept]
    write_tsv(clinical.groupby(["histology", "keep_pdac"]).size().rename("paired_cases").reset_index(),
              outdir / "cptac_histology_aggregate.tsv")
    raw.append(cptac)
    apgi_all = load_apgi(inputs / "apgi" / "protein_matrix.csv", inputs / "hgnc_complete_set.txt")
    write_tsv(pd.DataFrame(apgi_all["eligibility"]), outdir / "apgi_strict_pair_eligibility.tsv")
    strict_ix = apgi_all["strict_patient_indices"]
    strict_rec = {**apgi_all, "name": "APGI", "Dm": apgi_all["Dm"][:, strict_ix],
                  "pair_names": [apgi_all["pair_names"][i] for i in strict_ix]}
    mixed_rec = {**apgi_all, "name": "APGI_mixed", "Dm": apgi_all["Dm"],
                 "pair_names": list(apgi_all["pair_names"])}
    raw.append(strict_rec)
    raw.append(mixed_rec)
    raw.append(load_ruijin(inputs / "node_ruijin" / "20221227_PDAC_PRO_exp.xlsx"))
    proteins, dedup_summary, dedup_details = {}, [], []
    for rec in raw:
        c = select_highest_mean(rec["name"], rec["raw_symbols"], rec["abundance_X"], rec["Dm"], rec["source_labels"])
        c["pair_names"] = rec["pair_names"]
        c["pairs"] = len(rec["pair_names"])
        c["last_Dm"] = c.pop("last_Dm")
        proteins[rec["name"]] = c
        dedup_summary.append({**c["summary"], "paired_patients": len(rec["pair_names"])})
        dedup_details.extend(c["details"])
        log("protein cohort loaded", rec["name"], f"pairs={c['pairs']}",
            f"source_rows={c['summary']['rows_with_symbol']}", f"unique_symbols={len(c['genes'])}",
            f"duplicate_symbols={c['summary']['duplicate_symbols']}")
    write_tsv(pd.DataFrame(dedup_summary), outdir / "protein_dedup_summary.tsv")
    write_tsv(pd.DataFrame(dedup_details, columns=["cohort", "gene", "source_rows", "selected_row_index0",
              "selected_row_label", "selected_mean_abundance", "last_row_index0", "last_row_label",
              "last_mean_abundance", "selected_is_last_row", "selection_rule"]),
              outdir / "protein_dedup_duplicate_detail.tsv")

    # Compute whole-protein-universe paired statistics and an alternate last-row
    # selection solely as an audit of the historical dictionary-overwrite path.
    cohort_order = ["Fudan", "CPTAC", "APGI", "Ruijin"]
    pstats, laststats = {}, {}
    per_gene_rows = []
    lastrow_audit = []
    for name in cohort_order:
        c = proteins[name]
        ps = paired_stats(c["Dm"])
        ls = paired_stats(c["last_Dm"])
        pstats[name], laststats[name] = ps, ls
        for i, gene in enumerate(c["genes"]):
            per_gene_rows.append(dict(cohort=name, gene=str(gene), mean_delta=ps["mean"][i], sd_delta=ps["sd"][i],
                se_delta=ps["se"][i], p_value=ps["p"][i], q_value=ps["q"][i], n_pairs=int(ps["n"][i]),
                zero_sd_nonzero=bool(ps["sd"][i] == 0 and ps["n"][i] >= MIN_PAIRS and
                                     np.isfinite(ps["mean"][i]) and ps["mean"][i] != 0),
                zero_sd_zero_mean=bool(ps["sd"][i] == 0 and ps["n"][i] >= MIN_PAIRS and
                                       np.isfinite(ps["mean"][i]) and ps["mean"][i] == 0),
                mean_abundance=c["abundance"][i], testable=bool(ps["n"][i] >= MIN_PAIRS and np.isfinite(ps["mean"][i]))))
        for setname, frozen in (("R_protocol", R), ("R_strict", R5)):
            pos = {g: i for i, g in enumerate(c["genes"])}
            lpos = {g: i for i, g in enumerate(c["genes"])}
            total = len(frozen)
            mapped = [g for g in frozen if g in pos]
            ix = np.array([pos[g] for g in mapped], int)
            lix = np.array([lpos[g] for g in mapped], int)
            dirs = np.array([frozen[g] for g in mapped], int)
            test = (ps["n"][ix] >= MIN_PAIRS) & np.isfinite(ps["mean"][ix])
            sup = test & (ps["q"][ix] < PRIMARY_Q) & (np.sign(ps["mean"][ix]) == dirs)
            ltest = (ls["n"][lix] >= MIN_PAIRS) & np.isfinite(ls["mean"][lix])
            lsup = ltest & (ls["q"][lix] < PRIMARY_Q) & (np.sign(ls["mean"][lix]) == dirs)
            lastrow_audit.append(dict(cohort=name, frozen_set=setname, frozen_total=total,
                mapped_highest_mean=len(mapped), testable_highest_mean=int(test.sum()), supported_highest_mean=int(sup.sum()),
                mapped_last_row=len(mapped), testable_last_row=int(ltest.sum()), supported_last_row=int(lsup.sum()),
                supported_delta_last_minus_highest=int(lsup.sum() - sup.sum()),
                note="last-row selection is an audit comparator; primary uses max mean abundance"))
    write_tsv(pd.DataFrame(per_gene_rows), outdir / "protein_gene_paired_stats.tsv")
    bh_family_rows = []
    for name in cohort_order:
        st = pstats[name]
        bh_family_rows.append(dict(cohort=name, unique_mapped_protein_genes=len(proteins[name]["genes"]),
            estimable_BH_tests=int(np.isfinite(st["p"]).sum()),
            excluded_p_na=int((~np.isfinite(st["p"])).sum()),
            zero_sd_nonzero=int(np.sum((st["n"] >= MIN_PAIRS) & (st["sd"] == 0) &
                                       np.isfinite(st["mean"]) & (st["mean"] != 0))),
            zero_sd_zero_mean=int(np.sum((st["n"] >= MIN_PAIRS) & (st["sd"] == 0) &
                                         np.isfinite(st["raw_mean"]) & (st["raw_mean"] == 0))),
            min_pairs=MIN_PAIRS, BH_family="all finite protein paired-test p-values in this cohort"))
    write_tsv(pd.DataFrame(bh_family_rows), outdir / "protein_bh_family_summary.tsv")
    mixed = proteins["APGI_mixed"]
    mixed_stats = paired_stats(mixed["Dm"])
    write_tsv(pd.DataFrame([dict(cohort="APGI_mixed", gene=str(g), mean_delta=mixed_stats["mean"][i],
        sd_delta=mixed_stats["sd"][i], se_delta=mixed_stats["se"][i], p_value=mixed_stats["p"][i],
        q_value=mixed_stats["q"][i], n_pairs=int(mixed_stats["n"][i]),
        zero_sd_nonzero=bool(mixed_stats["sd"][i] == 0 and mixed_stats["n"][i] >= MIN_PAIRS and
                              np.isfinite(mixed_stats["mean"][i]) and mixed_stats["mean"][i] != 0),
        mean_abundance=mixed["abundance"][i], testable=bool(mixed_stats["n"][i] >= MIN_PAIRS and
                              np.isfinite(mixed_stats["mean"][i]))) for i,g in enumerate(mixed["genes"])]),
        outdir / "apgi_mixed_protein_gene_paired_stats.tsv")
    write_tsv(pd.DataFrame(lastrow_audit), outdir / "dedup_lastrow_sensitivity.tsv")
    with (outdir / "protein_patient_delta_long.tsv").open("w", encoding="utf-8", newline="") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(["cohort", "gene", "patient_id", "delta", "selected_source_row_index", "source_label"])
        for name in cohort_order + ["APGI_mixed"]:
            c = proteins[name]
            for gi, gene in enumerate(c["genes"]):
                for pi, patient in enumerate(c["pair_names"]):
                    w.writerow([name, str(gene), patient, c["Dm"][gi, pi],
                                int(c["selected_source_row_index"][gi]), c["source_labels"][gi]])

    # Build one authoritative long ledger over the shared A/B RNA universe.
    ledger_rows = []
    rank_arrays, score_arrays, dir_arrays = {}, {}, {}
    for rule in RULES:
        info = orders[rule]
        rnk = np.empty(len(master), dtype=int)
        rnk[info["order"]] = np.arange(1, len(master) + 1)
        rank_arrays[rule], score_arrays[rule], dir_arrays[rule] = rnk, info["score"], info["direction"]
    Rdir = np.array([R.get(str(g), 0) for g in master], int)
    Rset = np.array([str(g) in R for g in master], bool)
    R5set = np.array([str(g) in R5 for g in master], bool)
    for name in cohort_order:
        c, ps = proteins[name], pstats[name]
        pos = {g: i for i, g in enumerate(c["genes"])}
        pi = np.array([pos.get(str(g), -1) for g in master], int)
        mapped = pi >= 0
        n = np.full(len(master), np.nan)
        pm = np.full(len(master), np.nan)
        pp = np.full(len(master), np.nan)
        pq = np.full(len(master), np.nan)
        pab = np.full(len(master), np.nan)
        n[mapped] = ps["n"][pi[mapped]]
        pm[mapped] = ps["mean"][pi[mapped]]
        pp[mapped] = ps["p"][pi[mapped]]
        pq[mapped] = ps["q"][pi[mapped]]
        pab[mapped] = c["abundance"][pi[mapped]]
        testable = mapped & (n >= MIN_PAIRS) & np.isfinite(pm)
        support_R = Rset & testable & (pq < PRIMARY_Q) & (np.sign(pm) == Rdir)
        support_R5 = R5set & testable & (pq < PRIMARY_Q) & (np.sign(pm) == Rdir)
        base = dict(cohort=name, gene=master.astype(str), d_A=dA, se_A=seA, q_A=qA, d_B=dB, se_B=seB, q_B=qB,
                    in_R=Rset, in_R_strict=R5set, frozen_R_direction=Rdir, protein_mapped=mapped,
                    protein_mean_abundance=pab, protein_delta=pm, protein_p=pp, protein_q=pq,
                    protein_n_pairs=n, protein_testable=testable, supported_R=support_R, supported_R_strict=support_R5)
        df = pd.DataFrame(base)
        for rule in RULES:
            df[f"score_{rule}"] = score_arrays[rule]
            df[f"direction_{rule}"] = dir_arrays[rule]
            df[f"rank_{rule}"] = rank_arrays[rule]
            df[f"supported_{rule}"] = testable & (pq < PRIMARY_Q) & (np.sign(pm) == dir_arrays[rule])
        ledger_rows.append(df)
    ledger = pd.concat(ledger_rows, ignore_index=True)
    write_tsv(ledger, outdir / "analysis_ledger.tsv")
    mapping = ledger[["cohort", "gene", "protein_mapped", "protein_n_pairs", "protein_testable",
                      "protein_delta", "protein_p", "protein_q"]].copy()
    mapping["reason"] = np.select([~mapping.protein_mapped, mapping.protein_n_pairs < MIN_PAIRS,
                                   ~mapping.protein_testable],
                                  ["not_mapped_to_release", "insufficient_complete_pairs", "non_estimable"],
                                  default="testable")
    write_tsv(mapping, outdir / "protein_mapping_testability.tsv")
    log("authoritative ledger written", outdir / "analysis_ledger.tsv", f"rows={len(ledger)}")

    # Frozen-set coverage/testable/supported trio, all derived from the ledger.
    support_rows = []
    for name in cohort_order:
        subc = ledger[ledger.cohort == name]
        for setname, flag, supported_col in (("R_protocol", "in_R", "supported_R"),
                                             ("R_strict", "in_R_strict", "supported_R_strict")):
            sub = subc[subc[flag]]
            total = len(sub)
            mapped = int(sub.protein_mapped.sum())
            testable = int(sub.protein_testable.sum())
            supported = int(sub[supported_col].sum())
            support_rows.append(dict(cohort=name, set=setname, total_frozen_genes=total, mapped=mapped,
                testable_ge10=testable, supported_q05=supported, coverage=mapped / total if total else np.nan,
                testable_rate=testable / total if total else np.nan,
                unconditional_supported_rate=supported / total if total else np.nan,
                conditional_supported_rate=supported / testable if testable else np.nan,
                patients=proteins[name]["pairs"], min_pairs=MIN_PAIRS))
    support_df = pd.DataFrame(support_rows)
    write_tsv(support_df, outdir / "four_cohort_support.tsv")
    # Explicit 51-pair mixed-histology sensitivity; it never enters the strict
    # four-cohort shared-testable universe above.
    mixed_pos = {g: i for i, g in enumerate(mixed["genes"])}
    mix_pi = np.array([mixed_pos.get(str(g), -1) for g in master], int)
    mix_mapped = mix_pi >= 0
    mix_mean = np.full(len(master), np.nan); mix_p = np.full(len(master), np.nan)
    mix_q = np.full(len(master), np.nan); mix_n = np.full(len(master), np.nan)
    mix_mean[mix_mapped] = mixed_stats["mean"][mix_pi[mix_mapped]]
    mix_p[mix_mapped] = mixed_stats["p"][mix_pi[mix_mapped]]
    mix_q[mix_mapped] = mixed_stats["q"][mix_pi[mix_mapped]]
    mix_n[mix_mapped] = mixed_stats["n"][mix_pi[mix_mapped]]
    mix_test = mix_mapped & (mix_n >= MIN_PAIRS) & np.isfinite(mix_mean)
    mix_support_rows = []
    for setname, members in (("R_protocol", Rset), ("R_strict", R5set)):
        denom = int(members.sum()); mapped_n = int((members & mix_mapped).sum())
        test_n = int((members & mix_test).sum())
        sup_n = int((members & mix_test & (mix_q < PRIMARY_Q) & (np.sign(mix_mean) == Rdir)).sum())
        mix_support_rows.append(dict(cohort="APGI_mixed", set=setname, total_frozen_genes=denom,
            mapped=mapped_n, testable_ge10=test_n, supported_q05=sup_n,
            coverage=mapped_n / denom if denom else np.nan, testable_rate=test_n / denom if denom else np.nan,
            unconditional_supported_rate=sup_n / denom if denom else np.nan,
            conditional_supported_rate=sup_n / test_n if test_n else np.nan,
            patients=mixed["pairs"], min_pairs=MIN_PAIRS,
            cohort_label="mixed histology sensitivity: 44 PDAC + 4 Mucinous + 3 Others"))
    mix_support_df = pd.DataFrame(mix_support_rows)
    write_tsv(mix_support_df, outdir / "apgi_mixed_four_cohort_support.tsv")
    mix_ledger = pd.DataFrame(dict(cohort="APGI_mixed", gene=master.astype(str), protein_mapped=mix_mapped,
        protein_mean_delta=mix_mean, protein_p=mix_p, protein_q=mix_q, protein_n_pairs=mix_n,
        protein_testable=mix_test, testability_reason=np.where(~mix_mapped, "not_mapped",
            np.where(mix_n < MIN_PAIRS, "fewer_than_10_complete_pairs",
                     np.where(~np.isfinite(mix_mean), "non_estimable", "testable")))))
    for rule in RULES:
        mix_ledger[f"score_{rule}"] = score_arrays[rule]
        mix_ledger[f"direction_{rule}"] = dir_arrays[rule]
        mix_ledger[f"rank_{rule}"] = rank_arrays[rule]
        mix_ledger[f"supported_{rule}"] = mix_test & (mix_q < PRIMARY_Q) & (np.sign(mix_mean) == dir_arrays[rule])
    mix_ledger["in_R"] = Rset; mix_ledger["in_R_strict"] = R5set; mix_ledger["R_direction"] = Rdir
    write_tsv(mix_ledger, outdir / "apgi_mixed_analysis_ledger.tsv")
    mix_grid_rows = []
    for k in KS:
        for rule in RULES:
            sel = rank_arrays[rule] <= k
            mapped_n = int((sel & mix_mapped).sum()); test_n = int((sel & mix_test).sum())
            sup_n = int((sel & mix_test & (mix_q < PRIMARY_Q) &
                         (np.sign(mix_mean) == dir_arrays[rule])).sum())
            mix_grid_rows.append(dict(cohort="APGI_mixed", k=k, rule=rule, candidates=int(sel.sum()),
                mapped=mapped_n, testable=test_n, supported=sup_n, supported_per100=100*sup_n/k,
                testable_per100=100*test_n/k, conditional_supported_rate=sup_n/test_n if test_n else np.nan,
                denominator="ITT top-k; mixed-histology APGI 51-pair sensitivity"))
    write_tsv(pd.DataFrame(mix_grid_rows), outdir / "apgi_mixed_fixed_budget_grid.tsv")
    mix_dec = []
    for baseline in ("single_A", "single_B", "equal_weight"):
        bi = rank_arrays[baseline] <= 500; ci = rank_arrays["continuous"] <= 500
        bs = int((bi & mix_test & (mix_q < PRIMARY_Q) &
                  (np.sign(mix_mean) == dir_arrays[baseline])).sum())
        cs = int((ci & mix_test & (mix_q < PRIMARY_Q) &
                  (np.sign(mix_mean) == dir_arrays["continuous"])).sum())
        bt = int((bi & mix_test).sum()); ct = int((ci & mix_test).sum())
        rate = bs / bt if bt else 0.0; gain = cs - bs; coverage_component = (ct - bt) * rate
        mix_dec.append(dict(cohort="APGI_mixed", baseline=baseline, k=500,
            baseline_supported=bs, continuous_supported=cs, point_support_gain= gain,
            point_gain_per100=100*gain/500, baseline_testable=bt, continuous_testable=ct,
            point_testable_gain=ct-bt, baseline_conditional_rate=rate,
            point_coverage_component=coverage_component, point_residual=gain-coverage_component,
            interpretation="mixed-histology APGI sensitivity; descriptive decomposition"))
    write_tsv(pd.DataFrame(mix_dec), outdir / "apgi_mixed_descriptive_decomposition.tsv")
    gene_rng = np.random.default_rng(20261121)
    gene_ci_rows = []
    for name in cohort_order:
        subc = ledger[ledger.cohort == name]
        for setname, flag, supported_col in (("R_protocol", "in_R", "supported_R"),
                                             ("R_strict", "in_R_strict", "supported_R_strict")):
            sub = subc[subc[flag] & subc.protein_testable]
            outcomes = sub[supported_col].to_numpy(bool)
            draws = []
            if len(outcomes):
                for _ in range(1000):
                    take = gene_rng.integers(0, len(outcomes), len(outcomes))
                    draws.append(float(outcomes[take].mean()))
            gene_ci_rows.append(dict(cohort=name, set=setname, n_testable=len(outcomes),
                supported=int(outcomes.sum()), rate=float(outcomes.mean()) if len(outcomes) else np.nan,
                lo95=float(np.percentile(draws, 2.5)) if draws else np.nan,
                hi95=float(np.percentile(draws, 97.5)) if draws else np.nan,
                bootstrap_unit="testable genes sampled with replacement", draws=1000 if draws else 0,
                seed=20261121))
    write_tsv(pd.DataFrame(gene_ci_rows), outdir / "support_gene_bootstrap_ci.tsv")

    # Fixed-budget grid, with ITT denominator K: absent/untestable candidates fail.
    grid_rows = []
    for name in cohort_order:
        subc = ledger[ledger.cohort == name]
        for k in KS:
            for rule in RULES:
                sel = subc[subc[f"rank_{rule}"] <= k]
                mapped = int(sel.protein_mapped.sum())
                testable = int(sel.protein_testable.sum())
                supported = int(sel[f"supported_{rule}"].sum())
                grid_rows.append(dict(cohort=name, k=k, rule=rule, candidates=len(sel), mapped=mapped,
                    testable=testable, supported=supported, supported_per100=100.0 * supported / k,
                    testable_per100=100.0 * testable / k,
                    conditional_supported_rate=supported / testable if testable else np.nan,
                    denominator="ITT top-k; unmeasured or <10 pairs counts as unsupported"))
    grid_df = pd.DataFrame(grid_rows)
    write_tsv(grid_df, outdir / "fixed_budget_grid.tsv")

    # Descriptive, non-causal coverage decomposition at frozen k=500.
    dec_rows = []
    for name in cohort_order:
        subc = ledger[ledger.cohort == name]
        a = subc[subc.rank_single_A <= 500]
        cont = subc[subc.rank_continuous <= 500]
        a_s = int(a.supported_single_A.sum())
        c_s = int(cont.supported_continuous.sum())
        a_t, c_t = int(a.protein_testable.sum()), int(cont.protein_testable.sum())
        a_rate = a_s / a_t if a_t else 0.0
        extra = c_t - a_t
        coverage = extra * a_rate
        total = c_s - a_s
        dec_rows.append(dict(cohort=name, k=500, single_A_supported=a_s, continuous_supported=c_s,
            total_supported_gain=total, total_gain_per100=100.0 * total / 500,
            single_A_testable=a_t, continuous_testable=c_t, extra_testable=extra,
            single_A_conditional_rate=a_rate, coverage_component_supported=coverage,
            coverage_component_per100=100.0 * coverage / 500,
            residual_supported=total - coverage, residual_per100=100.0 * (total - coverage) / 500,
            interpretation="descriptive decomposition; not causal"))
    decomp_df = pd.DataFrame(dec_rows)
    write_tsv(decomp_df, outdir / "descriptive_decomposition.tsv")

    # A common all-cohort-testable candidate universe is a post hoc conditional
    # comparator. It supplements, never replaces, the unrestricted ITT grid.
    common = np.ones(len(master), bool)
    for name in cohort_order:
        sub = ledger[ledger.cohort == name].sort_values("gene")
        common &= sub.protein_testable.to_numpy(bool)
    common_genes = master[common]
    common_rank_rows = []
    common_support_rows = []
    for rule in RULES:
        order_common = orders[rule]["order"][common[orders[rule]["order"]]]
        selected = order_common[:500]
        for rank, i in enumerate(selected, 1):
            common_rank_rows.append(dict(rule=rule, rank=rank, gene=str(master[i]),
                                         score=score_arrays[rule][i], direction=int(np.sign(dir_arrays[rule][i])),
                                         comparator_role="conditional_posthoc_sensitivity"))
        for name in cohort_order:
            sub = ledger[(ledger.cohort == name) & (ledger.gene.isin(master[selected].astype(str)))]
            supported = int(sub[f"supported_{rule}"].sum())
            ksel = len(selected)
            common_support_rows.append(dict(cohort=name, rule=rule, common_testable_universe_n=len(common_genes),
                candidates=ksel, supported=supported, supported_per100=100.0 * supported / ksel if ksel else np.nan,
                conditional_supported_rate=supported / ksel if ksel else np.nan,
                comparator_role="conditional_posthoc_sensitivity", primary_ITT_substitute=False,
                shared_candidate_set_across_cohorts=True))
    write_tsv(pd.DataFrame(common_rank_rows), outdir / "common_testable_top500_lists.tsv")
    write_tsv(pd.DataFrame(common_support_rows), outdir / "common_testable_sensitivity.tsv")

    # Uniform protein-significance threshold sensitivity (all four cohorts).
    threshold_rows = []
    for name in cohort_order:
        subc = ledger[ledger.cohort == name]
        for label, cutoff, use_p in (("q<0.01", 0.01, False), ("q<0.05_primary", 0.05, False),
                                     ("q<0.10", 0.10, False), ("p<0.05_nominal", 0.05, True)):
            for setname, flag in (("R_protocol", "in_R"), ("R_strict", "in_R_strict")):
                s = subc[subc[flag]]
                sig = s.protein_p < cutoff if use_p else s.protein_q < cutoff
                sup = s.protein_testable & sig & (np.sign(s.protein_delta) == s.frozen_R_direction)
                threshold_rows.append(dict(scope="frozen_set", cohort=name, set=setname,
                    threshold=label, total=len(s), mapped=int(s.protein_mapped.sum()),
                    testable=int(s.protein_testable.sum()), supported=int(sup.sum()),
                    unconditional_rate=float(sup.sum() / len(s)) if len(s) else np.nan,
                    conditional_rate=float(sup.sum() / s.protein_testable.sum()) if s.protein_testable.sum() else np.nan,
                    consistent_across_all_four_cohorts=True))
            for rule in RULES:
                s = subc[subc[f"rank_{rule}"] <= 500]
                sig = s.protein_p < cutoff if use_p else s.protein_q < cutoff
                sup = s.protein_testable & sig & (np.sign(s.protein_delta) == s[f"direction_{rule}"])
                threshold_rows.append(dict(scope="unrestricted_top500", cohort=name, set=rule,
                    threshold=label, total=len(s), mapped=int(s.protein_mapped.sum()),
                    testable=int(s.protein_testable.sum()), supported=int(sup.sum()),
                    unconditional_rate=float(sup.sum() / len(s)) if len(s) else np.nan,
                    conditional_rate=float(sup.sum() / s.protein_testable.sum()) if s.protein_testable.sum() else np.nan,
                    consistent_across_all_four_cohorts=True))
    write_tsv(pd.DataFrame(threshold_rows), outdir / "protein_threshold_sensitivity.tsv")

    # Pair-count and magnitude filters are applied after each full-cohort BH
    # family has been calculated. Magnitude rows are exploratory only.
    sensitivity_rows = []
    for name in cohort_order:
        c = proteins[name]
        pos = {g: i for i, g in enumerate(c["genes"])}
        take = np.array([pos.get(str(g), -1) for g in master], int)
        mapped = take >= 0
        for min_pairs in (5, 10, 15, 20):
            st = paired_stats(c["Dm"], min_pairs=min_pairs)
            means = np.full(len(master), np.nan); qs = np.full(len(master), np.nan)
            ns = np.zeros(len(master), int)
            means[mapped] = st["mean"][take[mapped]]; qs[mapped] = st["q"][take[mapped]]
            ns[mapped] = st["n"][take[mapped]]
            testable = mapped & (ns >= min_pairs) & np.isfinite(means)
            for scope, selector, direction, setname in (
                [("frozen_set", Rset, Rdir, "R_protocol"), ("frozen_set", R5set, Rdir, "R_strict")] +
                [("unrestricted_top500", rank_arrays[r] <= 500, dir_arrays[r], r) for r in RULES]):
                denom = int(selector.sum())
                mapped_n = int((selector & mapped).sum()); test_n = int((selector & testable).sum())
                for floor in (0.0, 0.25, 0.5):
                    hit = selector & testable & (qs < PRIMARY_Q) & (np.sign(means) == direction)
                    if floor > 0:
                        hit &= np.abs(means) >= floor
                    supported = int(hit.sum())
                    sensitivity_rows.append(dict(cohort=name, scope=scope, set=setname, min_pairs=min_pairs,
                        effect_floor=floor, total=denom, mapped=mapped_n, testable=test_n, supported=supported,
                        supported_per100=100*supported/denom if denom else np.nan,
                        unconditional_rate=supported/denom if denom else np.nan,
                        conditional_rate=supported/test_n if test_n else np.nan,
                        BH_q_threshold=PRIMARY_Q, effect_filter_label=("none" if floor == 0 else "exploratory_BH_gated_effect_filter_no_threshold_FDR"),
                        effect_scale="input analyzed scale; paired mean delta"))
    sensitivity_df = pd.DataFrame(sensitivity_rows)
    write_tsv(sensitivity_df, outdir / "protein_sensitivity_grid.tsv")
    write_tsv(sensitivity_df[sensitivity_df.effect_floor == 0], outdir / "protein_sensitivity_min_pairs.tsv")
    write_tsv(sensitivity_df[(sensitivity_df.min_pairs == MIN_PAIRS) & (sensitivity_df.effect_floor > 0)],
              outdir / "protein_sensitivity_effect_threshold.tsv")

    # Historical comparisons are explicit references from the frozen 2026-11-21
    # report, never used to define or overwrite the current output.
    historical_support = {
        "Fudan": {"R_protocol_testable": 2235, "R_protocol_supported": 1533, "R_strict_testable": 553, "R_strict_supported": 456},
        "CPTAC": {"R_protocol_testable": 3456, "R_protocol_supported": 2097, "R_strict_testable": 902, "R_strict_supported": 752},
        "APGI": {"R_protocol_testable": 1464, "R_protocol_supported": 678, "R_strict_testable": 391, "R_strict_supported": 258},
        "Ruijin": {"R_protocol_testable": 1591, "R_protocol_supported": 1060, "R_strict_testable": 343, "R_strict_supported": 319},
    }
    historical_grid = {
        "Fudan": {"single_A": 192, "single_B": 206, "equal_weight": 231, "continuous": 229, "double_sig": 229},
        "CPTAC": {"single_A": 340, "single_B": 348, "equal_weight": 374, "continuous": 375, "double_sig": 378},
        "APGI": {"single_A": 109, "single_B": 105, "equal_weight": 126, "continuous": 129, "double_sig": 130},
        "Ruijin": {"single_A": 117, "single_B": 132, "equal_weight": 145, "continuous": 148, "double_sig": 141},
    }
    comparisons = []
    for row in support_rows:
        h = historical_support[row["cohort"]]
        if row["set"] == "R_protocol":
            metric_pairs = (("testable_R", "R_protocol_testable", "testable_ge10"),
                            ("supported_R", "R_protocol_supported", "supported_q05"))
        else:
            metric_pairs = (("testable_R_strict", "R_strict_testable", "testable_ge10"),
                            ("supported_R_strict", "R_strict_supported", "supported_q05"))
        for metric, histkey, curkey in metric_pairs:
            hv, cv = h[histkey], row[curkey]
            comparisons.append(dict(cohort=row["cohort"], analysis=row["set"], metric=metric,
                                    historical_20261121=hv, current_highest_mean=cv, delta=cv - hv,
                                    note="current values are recomputed after highest-mean duplicate collapse"))
    for row in grid_rows:
        if row["k"] == 500:
            hv = historical_grid[row["cohort"]][row["rule"]]
            comparisons.append(dict(cohort=row["cohort"], analysis="top500", metric=f"supported_{row['rule']}",
                                    historical_20261121=hv, current_highest_mean=row["supported"],
                                    delta=row["supported"] - hv,
                                    note="unrestricted ITT top-500; duplicate handling and deterministic ties audited"))
    write_tsv(pd.DataFrame(comparisons), outdir / "historical_comparison.tsv")

    # Paired-patient bootstrap: a shared patient resample is used for all genes
    # within each cohort and every ranking rule in that replicate.
    rng = np.random.default_rng(SEED)
    contrasts = {name: defaultdict(list) for name in cohort_order}
    common_draws = {name: defaultdict(list) for name in cohort_order}
    contrast_pairs = (("continuous_minus_single_A", "continuous", "single_A"),
                      ("continuous_minus_single_B", "continuous", "single_B"),
                      ("continuous_minus_equal_weight", "continuous", "equal_weight"),
                      ("double_sig_minus_continuous", "double_sig", "continuous"),
                      ("equal_weight_minus_single_A", "equal_weight", "single_A"),
                      ("equal_weight_minus_single_B", "equal_weight", "single_B"),
                      ("single_B_minus_single_A", "single_B", "single_A"))
    common_contrasts = {name: defaultdict(list) for name in cohort_order}
    main_draw_rows, common_draw_rows = [], []
    fixed_rate_draw_rows, fixed_rank_draw_rows, decomposition_draw_rows, zero_sd_boot_rows = [], [], [], []
    fixed_rank_contrasts = (("continuous_minus_single_A", "continuous", "single_A"),
                            ("continuous_minus_single_B", "continuous", "single_B"),
                            ("continuous_minus_equal_weight", "continuous", "equal_weight"),
                            ("double_sig_minus_continuous", "double_sig", "continuous"))
    fixed_rank_selected = {rule: orders[rule]["order"][:500] for rule in RULES}
    # map protein gene names once to reduce bootstrap work at selection time.
    ppos = {name: {g: i for i, g in enumerate(proteins[name]["genes"])} for name in cohort_order}
    master_protein_indices = {name: np.array([ppos[name].get(str(g), -1) for g in master], int) for name in cohort_order}
    for b in range(BOOTSTRAPS):
        ia_draw = rng.integers(0, A["Dm"].shape[1], A["Dm"].shape[1])
        ib_draw = rng.integers(0, B["Dm"].shape[1], B["Dm"].shape[1])
        sa = paired_stats(A["Dm"][:, ia_draw])
        sb = paired_stats(B["Dm"][:, ib_draw])
        ba, bb = sa["mean"][selA], sb["mean"][selB]
        bsea, bseb = sa["se"][selA], sb["se"][selB]
        bqa, bqb = sa["q"][selA], sb["q"][selB]
        b_orders = scores_and_orders(master, ba, bsea, bqa, bb, bseb, bqb)
        selected = {rule: b_orders[rule]["order"][:500] for rule in RULES}
        boot_prot = {}
        common_b = np.ones(len(master), bool)
        for name in cohort_order:
            c = proteins[name]
            idx = rng.integers(0, c["Dm"].shape[1], c["Dm"].shape[1])
            st = bootstrap_stats(c, idx)
            posidx = master_protein_indices[name]
            has = posidx >= 0
            pm = np.full(len(master), np.nan)
            pq = np.full(len(master), np.nan)
            pn = np.zeros(len(master), int)
            pm[has], pq[has], pn[has] = st["mean"][posidx[has]], st["q"][posidx[has]], st["n"][posidx[has]]
            test = has & (pn >= MIN_PAIRS) & np.isfinite(pm)
            boot_prot[name] = (pm, pq, pn, test)
            zero_sd_boot_rows.append(dict(replicate=b + 1, cohort=name,
                estimable_BH_tests=int(np.isfinite(st["p"]).sum()),
                zero_sd_nonzero=int(np.sum((st["n"] >= MIN_PAIRS) & (st["sd"] == 0) &
                                           np.isfinite(st["mean"]) & (st["mean"] != 0))),
                zero_sd_zero_mean=int(np.sum((st["n"] >= MIN_PAIRS) & (st["sd"] == 0) &
                                             np.isfinite(st["raw_mean"]) & (st["raw_mean"] == 0))),
                bootstrap_unit="protein paired patients; same replicate draw used for every gene", seed=SEED))
            common_b &= test
        for name in cohort_order:
            pm, pq, pn, test = boot_prot[name]
            c = proteins[name]
            posidx = master_protein_indices[name]
            has = posidx >= 0
            for setname, members in (("R_protocol", Rset), ("R_strict", R5set)):
                denom = int(members.sum())
                mapped_n = int((members & has).sum())
                test_n = int((members & test).sum())
                support_n = int((members & test & (pq < PRIMARY_Q) & (np.sign(pm) == Rdir)).sum())
                vals = {"mapped_per_total": mapped_n/denom if denom else np.nan,
                        "testable_per_total": test_n/denom if denom else np.nan,
                        "supported_per_total": support_n/denom if denom else np.nan,
                        "supported_per_testable": support_n/test_n if test_n else np.nan}
                for rate_name, rate in vals.items():
                    fixed_rate_draw_rows.append(dict(replicate=b+1, cohort=name, frozen_set=setname,
                        rate_name=rate_name, rate=rate, mapped=mapped_n, testable=test_n,
                        supported=support_n, denominator=denom, min_pairs=MIN_PAIRS,
                        protein_q_threshold=PRIMARY_Q, seed=SEED))
            supports = {}
            for rule in RULES:
                idx = selected[rule]
                good = test[idx] & (pq[idx] < PRIMARY_Q) & (np.sign(pm[idx]) == b_orders[rule]["direction"][idx])
                supports[rule] = int(good.sum())
            for label, left, right in contrast_pairs:
                val = 100.0 * (supports[left] - supports[right]) / 500.0
                contrasts[name][label].append(val)
                main_draw_rows.append(dict(replicate=b + 1, cohort=name, contrast=label,
                                           delta_supported_per100=val, paired_patient_resampling=True,
                                           seed=SEED))
            fixed_counts, fixed_tests = {}, {}
            for rule in RULES:
                idx = fixed_rank_selected[rule]
                fixed_tests[rule] = int(test[idx].sum())
                good = test[idx] & (pq[idx] < PRIMARY_Q) & (np.sign(pm[idx]) == orders[rule]["direction"][idx])
                fixed_counts[rule] = int(good.sum())
            for label, left, right in fixed_rank_contrasts:
                delta = 100.0 * (fixed_counts[left] - fixed_counts[right]) / 500.0
                fixed_rank_draw_rows.append(dict(replicate=b+1, cohort=name, contrast=label,
                    delta_supported_per100=delta, continuous_supported=fixed_counts["continuous"],
                    comparator_supported=fixed_counts[right], bootstrap_unit="protein paired patients only; RNA ranks fixed",
                    seed=SEED))
            for baseline in ("single_A", "single_B", "equal_weight"):
                base_test, cont_test = fixed_tests[baseline], fixed_tests["continuous"]
                base_sup, cont_sup = fixed_counts[baseline], fixed_counts["continuous"]
                base_rate = base_sup / base_test if base_test else 0.0
                gain = cont_sup - base_sup
                coverage_component = (cont_test - base_test) * base_rate
                decomposition_draw_rows.append(dict(replicate=b+1, cohort=name, baseline=baseline,
                    testable_gain_per100=100*(cont_test-base_test)/500,
                    support_gain_per100=100*gain/500,
                    coverage_component_per100=100*coverage_component/500,
                    residual_per100=100*(gain-coverage_component)/500,
                    baseline_testable=base_test, continuous_testable=cont_test,
                    baseline_supported=base_sup, continuous_supported=cont_sup,
                    bootstrap_unit="shared protein paired-patient draws; fixed RNA ranks", seed=SEED))
            if common_b.sum() >= 500:
                for rule in RULES:
                    order_common = b_orders[rule]["order"][common_b[b_orders[rule]["order"]]][:500]
                    candidates = len(order_common)
                    good = (pq[order_common] < PRIMARY_Q) & (np.sign(pm[order_common]) == b_orders[rule]["direction"][order_common])
                    n_sup = int(good.sum())
                    rate = 100.0 * n_sup / candidates if candidates else np.nan
                    common_draws[name][rule].append(rate)
                    common_draw_rows.append(dict(replicate=b + 1, cohort=name, rule=rule,
                        common_testable_universe_n=int(common_b.sum()), candidates=candidates,
                        supported=n_sup, supported_per100=rate,
                        comparator_role="conditional_posthoc_sensitivity",
                        same_replicate_shared_across_methods_and_cohorts=True, seed=SEED))
                if common_draws[name]:
                    cb_supports = {rule: common_draws[name][rule][-1] for rule in RULES}
                    for label, left, right in contrast_pairs:
                        common_contrasts[name][label].append(cb_supports[left] - cb_supports[right])
            if b == 0:
                log("bootstrap started; first replicate completed")
        if (b + 1) % 100 == 0:
            log("bootstrap progress", f"{b + 1}/{BOOTSTRAPS}")

    write_tsv(pd.DataFrame(main_draw_rows), outdir / "paired_patient_bootstrap_draws.tsv")
    write_tsv(pd.DataFrame(common_draw_rows), outdir / "common_testable_bootstrap_draws.tsv")
    boot_rows = []
    for name in cohort_order:
        for label, vals in contrasts[name].items():
            v = np.asarray(vals, float)
            boot_rows.append(dict(cohort=name, contrast=label, draws=len(v), mean_delta_per100=float(v.mean()),
                median_delta_per100=float(np.median(v)), lo95=float(np.percentile(v, 2.5)),
                hi95=float(np.percentile(v, 97.5)), probability_delta_gt0=float(np.mean(v > 0)),
                bootstrap_unit="paired patient deltas resampled within RNA and each protein cohort",
                seed=SEED, primary_scope="unrestricted top-500 ITT"))
    write_tsv(pd.DataFrame(boot_rows), outdir / "paired_patient_bootstrap_summary.tsv")
    common_boot_rows = []
    for name in cohort_order:
        for rule, vals in common_draws[name].items():
            v = np.asarray(vals, float)
            common_boot_rows.append(dict(cohort=name, rule=rule, draws=len(v), mean_supported_per100=float(v.mean()),
                lo95=float(np.percentile(v, 2.5)), hi95=float(np.percentile(v, 97.5)),
                comparator_role="conditional_posthoc_sensitivity", seed=SEED))
        for label, vals in common_contrasts[name].items():
            v = np.asarray(vals, float)
            common_boot_rows.append(dict(cohort=name, rule=label, draws=len(v), mean_supported_per100=float(v.mean()),
                lo95=float(np.percentile(v, 2.5)), hi95=float(np.percentile(v, 97.5)),
                comparator_role="conditional_posthoc_paired_contrast", seed=SEED))
    write_tsv(pd.DataFrame(common_boot_rows), outdir / "common_testable_bootstrap_summary.tsv")

    def summarize(rows, group_keys, point, value_key):
        records = []
        frame = pd.DataFrame(rows)
        for keys, block in frame.groupby(group_keys, sort=False):
            values = block[value_key].dropna().to_numpy(float)
            lo, hi = np.percentile(values, (2.5, 97.5))
            records.append(dict(zip(group_keys, keys), point=point[keys], mean=float(values.mean()),
                                median=float(np.median(values)), lo95=lo, hi95=hi,
                                draws=len(values), total_draws=BOOTSTRAPS, seed=SEED))
        return pd.DataFrame(records)

    fixed_points = {}
    metric_names = {"mapped_per_total": "coverage", "testable_per_total": "testable_rate",
                    "supported_per_total": "unconditional_supported_rate",
                    "supported_per_testable": "conditional_supported_rate"}
    for row in support_df.itertuples(index=False):
        for metric, col in metric_names.items():
            fixed_points[row.cohort, row.set, metric] = getattr(row, col)
    fixed_set_summary = summarize(fixed_rate_draw_rows, ["cohort", "frozen_set", "rate_name"], fixed_points, "rate")
    fixed_set_summary = fixed_set_summary.rename(columns={"frozen_set": "set"})
    fixed_set_summary["estimand"] = "paired protein patients; RNA set and direction fixed"
    write_tsv(fixed_set_summary, outdir / "protein_patient_bootstrap_support.tsv")

    observed = grid_df[grid_df.k == 500].set_index(["cohort", "rule"])
    points = {(name, label): float(observed.loc[(name, left), "supported_per100"] -
                                  observed.loc[(name, right), "supported_per100"])
              for name in cohort_order for label, left, right in fixed_rank_contrasts}
    fixed_ranking = summarize(fixed_rank_draw_rows, ["cohort", "contrast"], points, "delta_supported_per100")
    fixed_ranking["estimand"] = "paired protein patients; RNA ranking and direction fixed"
    write_tsv(fixed_ranking, outdir / "fixed_rank_protein_bootstrap_summary.tsv")

    component_points = {}
    for name in cohort_order:
        continuous = observed.loc[(name, "continuous")]
        for baseline in ("single_A", "single_B", "equal_weight"):
            other = observed.loc[(name, baseline)]
            testable_gain = (continuous.testable - other.testable) / 5
            support_gain = (continuous.supported - other.supported) / 5
            coverage_component = testable_gain * other.conditional_supported_rate
            values = [testable_gain, support_gain, coverage_component, support_gain - coverage_component]
            metrics = ["testable_gain_per100", "support_gain_per100", "coverage_component_per100", "residual_per100"]
            for metric, value in zip(metrics, values):
                component_points[name, baseline, metric] = value
    component_draw = pd.DataFrame(decomposition_draw_rows).melt(
        id_vars=["cohort", "baseline", "replicate"], value_vars=metrics, var_name="component", value_name="delta")
    dec_summary = summarize(component_draw.to_dict("records"), ["cohort", "baseline", "component"], component_points, "delta")
    dec_summary["estimand"] = "paired protein patients; fixed rankings and observed baseline rate identity within each draw"
    write_tsv(dec_summary, outdir / "decomposition_bootstrap_summary.tsv")

    zero_sd = pd.DataFrame(zero_sd_boot_rows)
    write_tsv(zero_sd, outdir / "bootstrap_degenerate_tests.tsv")

    # Record cohort sizes and output hashes for reproducibility.
    cohort_meta = {name: dict(pairs=proteins[name]["pairs"], proteins=len(proteins[name]["genes"]),
        duplicate_symbols=dedup_summary[cohort_order.index(name)]["duplicate_symbols"],
        rows_removed=dedup_summary[cohort_order.index(name)]["rows_removed_by_dedup"])
        for name in cohort_order}
    metadata = dict(started_utc=datetime.fromtimestamp(t0, timezone.utc).isoformat(),
        completed_utc=datetime.now(timezone.utc).isoformat(), elapsed_seconds=round(time.time() - t0, 2),
        python=sys.version, numpy=np.__version__, pandas=pd.__version__, scipy=stats.__version__ if hasattr(stats, "__version__") else "see package metadata",
        seed=SEED, bootstrap_draws=BOOTSTRAPS, gene_support_ci_draws=1000,
        gene_support_ci_seed=20261121, min_pairs=MIN_PAIRS, protein_q_threshold=PRIMARY_Q,
        frozen_protocol_R_count=len(R), frozen_strict_R_count=len(R5), shared_RNA_gene_count=len(master),
        RNA_A_pairs=A["pairs"], RNA_B_pairs=B["pairs"], RNA_A_genes=len(A["genes"]), RNA_B_genes=len(B["genes"]),
        protein_cohorts=cohort_meta,
        ranking_rules=list(RULES), k_grid=list(KS), duplicate_rule="highest mean transformed/observed protein abundance; ties -> earliest source row",
        primary_support_direction="frozen R_protocol dir = sign(d_A + d_B)",
        bootstrap_resampling="paired delta columns resampled together within each cohort; one patient index draw applied to every gene in that cohort",
        common_testable_sensitivity="post hoc conditional comparator; not a substitute for unrestricted ITT",
        cptac_histology_mapping_sha256=sha256(args.cptac_histology),
        inputs=manifest)
    (outdir / "run_metadata.json").write_text(json.dumps(metadata, indent=2, default=str), encoding="utf-8")
    out_manifest = []
    for p in sorted(outdir.rglob("*")):
        if p.is_file() and p.name != "output_manifest.tsv":
            out_manifest.append(dict(file=str(p.relative_to(outdir)), bytes=p.stat().st_size, sha256=sha256(p)))
    write_tsv(pd.DataFrame(out_manifest), outdir / "output_manifest.tsv")
    log("DONE", f"elapsed_s={round(time.time() - t0, 1)}", f"output={outdir}",
        f"manifest={outdir / 'output_manifest.tsv'}")


if __name__ == "__main__":
    main()
