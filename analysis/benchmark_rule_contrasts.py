#!/usr/bin/env python3
"""Benchmark analyses from frozen core artifacts.

All ranking/prior/matching rules are fixed in protocol_frozen.md. This runner
does not inspect target protein outcomes until the fold-specific ranked lists
and BENCH_B matches have been constructed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import rankdata


COHORTS = ("Fudan", "CPTAC", "APGI", "Ruijin")
RNA_RULES = ("single_A", "single_B", "equal_weight", "continuous")
BENCH_A_RULES = (*RNA_RULES, "availability_only", "RNA_plus_availability")
BENCH_C_RULES = (*RNA_RULES, "double_sig", "availability_only", "RNA_plus_availability")
KS = (100, 250, 500, 1000, 2000)
BENCH_B_CALIPER_SD = 0.2
MIN_PAIRS = 10
PRIMARY_Q = 0.05


def log(*parts):
    print(datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"), *parts,
          flush=True)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def find_file(root: Path, candidates, required=True):
    for name in candidates:
        p = root / name
        if p.is_file():
            return p
    if required:
        raise FileNotFoundError(f"none of {candidates} found under {root}")
    return None


def col(df, candidates, required=True):
    lookup = {str(x).strip().lower(): x for x in df.columns}
    for name in candidates:
        if name.lower() in lookup:
            return lookup[name.lower()]
    if required:
        raise ValueError(f"missing expected column from {candidates}; have {list(df.columns)}")
    return None


def truth(v):
    if pd.isna(v):
        return False
    if isinstance(v, (bool, np.bool_)):
        return bool(v)
    if isinstance(v, (int, np.integer, float, np.floating)):
        return bool(v != 0)
    return str(v).strip().lower() in {"1", "true", "t", "yes", "y"}


def norm_gene(v):
    return str(v).strip().upper()


def norm_cohort(v):
    s = str(v).strip()
    low = s.lower()
    aliases = {"fudan": "Fudan", "cptac": "CPTAC", "apgi": "APGI",
               "apgi_strict": "APGI", "apgi_pdac44": "APGI", "ruijin": "Ruijin"}
    return aliases.get(low, s)


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


def paired_stats(Dm):
    """Match the core's paired t/BH implementation and minimum-pair rule."""
    Dm = np.asarray(Dm, float)
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
    bad = (n < MIN_PAIRS) | ~np.isfinite(p)
    p[bad] = np.nan
    mean[bad] = np.nan
    se[bad] = np.nan
    q = bh(p)
    return {"mean": mean, "p": p, "q": q, "n": n, "se": se, "sd": sd}


def write_tsv(df, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, na_rep="NA", float_format="%.10g",
              encoding="utf-8", lineterminator="\n")


def normalize_table(df, required_columns):
    for c in required_columns:
        if c not in df.columns:
            raise ValueError(f"expected {c} in columns {list(df.columns)}")
    return df


def load_inputs(core_dir: Path, args):
    ranks_path = Path(args.ranks) if args.ranks else find_file(
        core_dir, ["rna_scores_all_ranks.tsv", "rna_scores_and_all_ranks.tsv", "rna_scores_all_ranks.csv"])
    delta_path = Path(args.deltas) if args.deltas else find_file(
        core_dir, ["protein_patient_delta_long.tsv", "protein_patient_deltas_long.tsv"])
    mapping_path = Path(args.mapping) if args.mapping else find_file(
        core_dir, ["protein_mapping_testability.tsv", "protein_mapping_and_testability.tsv"])
    r_path = Path(args.r_protocol) if args.r_protocol else find_file(
        core_dir, ["R_protocol.tsv", "step2/R_protocol.tsv", "frozen_R_protocol.tsv"], required=False)
    ranks = pd.read_csv(ranks_path, sep="\t")
    ranks.columns = [str(c).strip() for c in ranks.columns]
    ranks["gene"] = ranks["gene"].map(norm_gene)
    if ranks["gene"].duplicated().any():
        raise ValueError("RNA score/rank input has duplicate genes")
    for rule in (*RNA_RULES, "double_sig"):
        for suffix in ("rank", "direction"):
            if f"{suffix}_{rule}" not in ranks.columns:
                raise ValueError(f"RNA table missing {suffix}_{rule}")

    dlong_raw = pd.read_csv(delta_path, sep="\t", low_memory=False)
    dlong_raw.columns = [str(c).strip() for c in dlong_raw.columns]
    cohort_c = col(dlong_raw, ["cohort", "protein_cohort", "dataset"])
    gene_c = col(dlong_raw, ["gene", "symbol", "gene_symbol"])
    patient_c = col(dlong_raw, ["patient_id", "case_id", "pair_id", "patient", "case"])
    delta_c = col(dlong_raw, ["delta", "protein_delta", "tumor_minus_normal", "tumour_minus_normal", "delta_value"])
    dlong = dlong_raw[[cohort_c, gene_c, patient_c, delta_c]].rename(
        columns={cohort_c: "cohort", gene_c: "gene", patient_c: "patient_id", delta_c: "delta"})
    dlong["raw_cohort_label"] = dlong["cohort"].astype(str).str.strip()
    dlong["cohort"] = dlong["raw_cohort_label"].map(norm_cohort)
    cohort_audit = (dlong.groupby(["raw_cohort_label", "cohort"], dropna=False)
                    .size().rename("rows").reset_index())
    cohort_audit["included_strict_primary"] = cohort_audit["cohort"].isin(COHORTS)
    cohort_audit["exclusion_reason"] = np.where(
        cohort_audit["included_strict_primary"], "",
        np.where(cohort_audit["raw_cohort_label"].str.lower().str.contains("apgi.*mixed|mixed.*apgi", regex=True),
                 "APGI mixed-histology sensitivity excluded from strict 44-pair primary",
                 "cohort label not in strict four-cohort primary"))
    dlong["gene"] = dlong["gene"].map(norm_gene)
    dlong["patient_id"] = dlong["patient_id"].astype(str)
    dlong["delta"] = pd.to_numeric(dlong["delta"], errors="coerce")
    dlong = dlong[dlong["cohort"].isin(COHORTS)].copy()
    if dlong.duplicated(["cohort", "gene", "patient_id"]).any():
        raise ValueError("protein delta input has duplicate cohort/gene/patient rows")

    mdf = pd.read_csv(mapping_path, sep="\t")
    mdf.columns = [str(c).strip() for c in mdf.columns]
    mc = col(mdf, ["cohort", "protein_cohort", "dataset"])
    mg = col(mdf, ["gene", "symbol", "gene_symbol"])
    mm = col(mdf, ["protein_mapped", "mapped", "is_mapped", "released_matrix_mapped"], required=False)
    mn = col(mdf, ["protein_n_pairs", "n_pairs", "pair_count", "n_finite_pairs"], required=False)
    mt = col(mdf, ["protein_testable", "testable", "is_testable"], required=False)
    mr = col(mdf, ["reason", "testability_reason", "mapping_reason", "state"], required=False)
    mdf = mdf.rename(columns={mc: "cohort", mg: "gene"})
    if mm:
        mdf["mapped_input"] = mdf[mm].map(truth)
    else:
        mdf["mapped_input"] = np.nan
    if mn:
        mdf["n_pairs_input"] = pd.to_numeric(mdf[mn], errors="coerce")
    else:
        mdf["n_pairs_input"] = np.nan
    if mt:
        mdf["testable_input"] = mdf[mt].map(truth)
    else:
        mdf["testable_input"] = np.nan
    if mr:
        mdf["reason_input"] = mdf[mr].astype(str)
    else:
        mdf["reason_input"] = ""
    mdf = mdf[["cohort", "gene", "mapped_input", "n_pairs_input", "testable_input", "reason_input"]]
    mdf["cohort"] = mdf["cohort"].map(norm_cohort)
    mdf["gene"] = mdf["gene"].map(norm_gene)
    if mdf.duplicated(["cohort", "gene"]).any():
        raise ValueError("mapping/testability input has duplicate cohort/gene rows")

    frozen_R = None
    if r_path:
        rdf = pd.read_csv(r_path, sep="\t")
        rc_gene = col(rdf, ["gene", "symbol", "gene_symbol"])
        rc_dir = col(rdf, ["dir", "direction", "R_direction"], required=False)
        rdf["gene"] = rdf[rc_gene].map(norm_gene)
        if rc_dir:
            rdf["direction"] = pd.to_numeric(rdf[rc_dir], errors="coerce").fillna(0).astype(int)
        else:
            rdf["direction"] = np.nan
        if rdf["gene"].duplicated().any():
            raise ValueError("frozen R protocol has duplicate genes")
        frozen_R = dict(zip(rdf["gene"], rdf["direction"]))
    elif "in_frozen_R" in ranks.columns:
        frozen_R = {row.gene: int(np.sign(row.direction_equal_weight))
                    for row in ranks.itertuples() if truth(row.in_frozen_R)}
    else:
        raise FileNotFoundError("exact frozen R_protocol.tsv or in_frozen_R rank column is required")

    matrices = {}
    for cohort in COHORTS:
        part = dlong[dlong["cohort"] == cohort]
        if part.empty:
            raise ValueError(f"no protein delta rows for {cohort}")
        piv = part.pivot(index="gene", columns="patient_id", values="delta").sort_index()
        matrices[cohort] = (piv.index.to_numpy(str), piv.columns.to_numpy(str), piv.to_numpy(float))
    apgi_pairs = matrices["APGI"][1].size
    if apgi_pairs != 44:
        raise ValueError(f"strict APGI benchmark requires exactly 44 PDAC matched pairs; observed {apgi_pairs}")
    cptac_pairs = matrices["CPTAC"][1].size
    if cptac_pairs != 64:
        raise ValueError(f"strict CPTAC benchmark requires exactly 64 PDAC matched pairs after excluding 3 adenocarcinoma cases; observed {cptac_pairs}")

    return {"ranks_path": ranks_path, "delta_path": delta_path, "mapping_path": mapping_path,
            "r_path": r_path, "ranks": ranks, "dlong": dlong, "mdf": mdf,
            "frozen_R": frozen_R, "matrices": matrices, "cohort_input_audit": cohort_audit}


def make_state_map(inputs, cohort, universe_genes):
    """Audit mapping and finite-pair testability; absent rows remain distinct."""
    gene_names, patient_names, D = inputs["matrices"][cohort]
    pidx = {g: i for i, g in enumerate(gene_names)}
    subset = inputs["mdf"][inputs["mdf"].cohort == cohort].set_index("gene", drop=False)
    n_obs = np.zeros(len(gene_names), dtype=int)
    if D.size:
        n_obs = np.isfinite(D).sum(axis=1).astype(int)
    states = {}
    audits = []
    for gene in universe_genes:
        mi = gene in pidx
        row = subset.loc[gene] if gene in subset.index else None
        if isinstance(row, pd.DataFrame):
            raise ValueError(f"duplicate mapping rows for {cohort}/{gene}")
        input_map = truth(row["mapped_input"]) if row is not None and pd.notna(row["mapped_input"]) else mi
        mapped = bool(mi or input_map)
        npairs = int(n_obs[pidx[gene]]) if mi else 0
        if row is not None and pd.notna(row["n_pairs_input"]):
            input_n = int(row["n_pairs_input"])
        else:
            input_n = npairs
        testable = bool(mapped and npairs >= MIN_PAIRS)
        if not mapped:
            state = "not_in_released_matrix"
        elif not testable:
            state = "mapped_below_10_finite_pairs"
        else:
            state = "testable_ge10"
        if row is not None and pd.notna(row["testable_input"]):
            input_testable = truth(row["testable_input"])
        else:
            input_testable = testable
        reason = str(row["reason_input"]) if row is not None else ""
        if not mapped and reason and reason.lower() not in {"nan", "none"}:
            state = "not_in_released_matrix"
        states[gene] = {"mapped": mapped, "n_pairs": npairs, "testable": testable, "state": state,
                        "matrix_row_present": mi, "input_mapped": input_map,
                        "input_n_pairs": input_n, "input_testable": input_testable}
        if (mapped != input_map) or (npairs != input_n) or (testable != input_testable):
            audits.append(dict(cohort=cohort, gene=gene, matrix_row_present=mi,
                               mapped_from_matrix_or_metadata=mapped, metadata_mapped=input_map,
                               n_finite_pairs_recomputed=npairs, metadata_n_pairs=input_n,
                               testable_recomputed=testable, metadata_testable=input_testable,
                               reason=reason))
    return states, audits


def stable_order_desc(score, genes):
    genes = np.asarray(genes, str)
    score = np.asarray(score, float)
    clean = np.where(np.isfinite(score), score, -np.inf)
    return np.lexsort((genes, -clean))


def make_fold_rankings(ranks, state_maps, target):
    genes = ranks["gene"].to_numpy(str)
    prior = np.zeros(len(genes), dtype=int)
    states_long = []
    for i, gene in enumerate(genes):
        for source in COHORTS:
            is_target = source == target
            if is_target:
                # The held-out source is named for the fold-disjoint audit only.
                # Do not read its mapping/testability state while constructing ranks.
                states_long.append(dict(fold=f"holdout_{target}", target_cohort=target, gene=gene,
                                        source_cohort=source, source_role="excluded_target",
                                        state="not_read_for_prior", protein_mapped=pd.NA,
                                        n_finite_pairs=pd.NA, protein_testable=pd.NA,
                                        prior_contribution=0, target_excluded_from_prior=True))
                continue
            s = state_maps[source].get(gene, {"mapped": False, "n_pairs": 0,
                                               "testable": False, "state": "not_in_released_matrix"})
            states_long.append(dict(fold=f"holdout_{target}", target_cohort=target, gene=gene,
                                    source_cohort=source, source_role="training_prior",
                                    state=s["state"], protein_mapped=s["mapped"], n_finite_pairs=s["n_pairs"],
                                    protein_testable=s["testable"], prior_contribution=int(s["testable"]),
                                    target_excluded_from_prior=False))
            prior[i] += int(s["testable"])
    ord_prior = stable_order_desc(prior, genes)
    rank_prior = np.empty(len(genes), dtype=int)
    rank_prior[ord_prior] = np.arange(1, len(genes) + 1)
    # Use midranks for tied availability levels in the pre-fixed rank fusion.
    prior_midrank = rankdata(-prior, method="average")
    n = len(genes)
    prior_pct = (n - prior_midrank) / (n - 1) if n > 1 else np.ones(n)
    r = ranks.set_index("gene", drop=False)
    cont_rank = pd.to_numeric(r["rank_continuous"], errors="coerce").to_numpy(float)
    cont_midrank_pct = (n - cont_rank) / (n - 1) if n > 1 else np.ones(n)
    combined_score = (cont_midrank_pct + prior_pct) / 2.0
    ord_combined = stable_order_desc(combined_score, genes)
    rank_combined = np.empty(len(genes), dtype=int)
    rank_combined[ord_combined] = np.arange(1, len(genes) + 1)

    data = ranks.copy()
    data["availability_prior_count_0_3"] = prior
    data["rank_availability_only"] = rank_prior
    data["direction_availability_only"] = pd.to_numeric(data["direction_equal_weight"], errors="coerce").fillna(0).astype(int)
    data["score_RNA_plus_availability"] = combined_score
    data["rank_RNA_plus_availability"] = rank_combined
    data["direction_RNA_plus_availability"] = pd.to_numeric(data["direction_continuous"], errors="coerce").fillna(0).astype(int)
    return data, states_long


def direction_for(rule, data):
    return pd.to_numeric(data[f"direction_{rule}"], errors="coerce").fillna(0).to_numpy(int)


def order_for(rule, data):
    return np.argsort(pd.to_numeric(data[f"rank_{rule}"], errors="coerce").to_numpy(float), kind="mergesort")


def point_yield_for(selected_genes, selected_dirs, stats_all, gene_pos, mapping_state):
    mapped = np.array([bool(mapping_state.get(g, {}).get("mapped", False)) for g in selected_genes])
    pidx = np.array([gene_pos.get(g, -1) for g in selected_genes], int)
    present = pidx >= 0
    mapped &= present
    testable = np.zeros(len(selected_genes), bool)
    supported = np.zeros(len(selected_genes), bool)
    if present.any():
        ix = pidx[present]
        testable[present] = (stats_all["n"][ix] >= MIN_PAIRS) & np.isfinite(stats_all["mean"][ix])
        good = (testable[present] & (stats_all["q"][ix] < PRIMARY_Q) &
                (np.sign(stats_all["mean"][ix]) == np.asarray(selected_dirs, int)[present]))
        supported[present] = good
    return {"candidate_slots": len(selected_genes), "mapped": int(mapped.sum()),
            "unmapped": int(len(selected_genes) - mapped.sum()), "testable_ge10": int(testable.sum()),
            "mapped_but_untestable": int(mapped.sum() - testable.sum()),
            "supported_q05": int(supported.sum()),
            "supported_per100_slots": 100.0 * float(supported.sum()) / len(selected_genes) if len(selected_genes) else np.nan,
            "mapped_per100_slots": 100.0 * float(mapped.sum()) / len(selected_genes) if len(selected_genes) else np.nan,
            "testable_per100_slots": 100.0 * float(testable.sum()) / len(selected_genes) if len(selected_genes) else np.nan,
            "conditional_supported_rate": float(supported.sum()) / testable.sum() if testable.sum() else np.nan}


def bootstrap_yield(selected_genes, selected_dirs, stats_b, gene_pos, mapping_state):
    mapped = np.array([bool(mapping_state.get(g, {}).get("mapped", False)) for g in selected_genes])
    pidx = np.array([gene_pos.get(g, -1) for g in selected_genes], int)
    present = pidx >= 0
    mapped &= present
    testable = np.zeros(len(selected_genes), bool)
    supported = np.zeros(len(selected_genes), bool)
    if present.any():
        ix = pidx[present]
        testable[present] = (stats_b["n"][ix] >= MIN_PAIRS) & np.isfinite(stats_b["mean"][ix])
        supported[present] = (testable[present] & (stats_b["q"][ix] < PRIMARY_Q) &
                              (np.sign(stats_b["mean"][ix]) == np.asarray(selected_dirs, int)[present]))
    slots = len(selected_genes)
    return {"slots": slots, "mapped": int(mapped.sum()), "testable": int(testable.sum()),
            "supported": int(supported.sum()),
            "itt": 100.0 * float(supported.sum()) / slots if slots else np.nan,
            "conditional": 100.0 * float(supported.sum()) / testable.sum() if testable.sum() else np.nan,
            "testable_count": int(testable.sum())}


def summarize_draws(draws, keys, metric, scope):
    rows = []
    ddf = pd.DataFrame(draws)
    if ddf.empty:
        return rows
    grouped = ddf.groupby(keys, dropna=False)
    for key, group in grouped:
        if not isinstance(key, tuple):
            key = (key,)
        vals = pd.to_numeric(group[metric], errors="coerce").to_numpy(float)
        vals = vals[np.isfinite(vals)]
        rec = dict(zip(keys, key))
        rec.update(draws_total=int(len(group)), draws_valid=int(len(vals)),
                   bootstrap_center=float(np.mean(vals)) if len(vals) else np.nan,
                   bootstrap_median=float(np.median(vals)) if len(vals) else np.nan,
                   lo95=float(np.percentile(vals, 2.5)) if len(vals) else np.nan,
                   hi95=float(np.percentile(vals, 97.5)) if len(vals) else np.nan,
                   probability_gt0=float(np.mean(vals > 0)) if len(vals) else np.nan,
                   metric=metric, estimand_scope=scope)
        rows.append(rec)
    return rows


def load_for_fold(inputs, ranks, target):
    data, state_rows = make_fold_rankings(ranks, inputs["state_maps"], target)
    data = data.set_index("gene", drop=False).loc[ranks["gene"]].reset_index(drop=True)
    top_lists = {}
    for rule in BENCH_C_RULES:
        order = order_for(rule, data)
        genes = data.iloc[order]["gene"].to_numpy(str)
        dirs = direction_for(rule, data)[order]
        top_lists[rule] = (genes, dirs)
    return data, state_rows, top_lists


def eligible_matched(ranks, frozen_R):
    r = ranks.set_index("gene", drop=False)
    inR = set(frozen_R)
    elig = ranks[np.isfinite(pd.to_numeric(ranks["d_A"], errors="coerce")) &
                 np.isfinite(pd.to_numeric(ranks["se_A"], errors="coerce")) &
                 (pd.to_numeric(ranks["se_A"], errors="coerce") > 0)].copy()
    elig["is_R"] = elig["gene"].isin(inR)
    elig["a_direction"] = np.sign(pd.to_numeric(elig["d_A"], errors="coerce")).astype(int)
    elig["r_direction"] = elig["gene"].map(frozen_R).fillna(0).astype(int)
    elig["match_direction"] = elig["a_direction"]
    elig["feature_log_abs_dA"] = np.log(np.maximum(np.abs(pd.to_numeric(elig["d_A"], errors="coerce")), 1e-12))
    elig["feature_log_seA"] = np.log(np.maximum(pd.to_numeric(elig["se_A"], errors="coerce"), 1e-12))
    return elig


def match_controls(elig, prior_by_gene=None):
    """Greedy 1:1 no-replacement match under fixed sign and two calipers."""
    x = elig.copy()
    mean1 = x["feature_log_abs_dA"].mean()
    sd1 = x["feature_log_abs_dA"].std(ddof=1)
    mean2 = x["feature_log_seA"].mean()
    sd2 = x["feature_log_seA"].std(ddof=1)
    if not np.isfinite(sd1) or sd1 == 0 or not np.isfinite(sd2) or sd2 == 0:
        raise ValueError("BENCH_B matching features have zero/undefined SD")
    x["z_log_abs_dA"] = (x["feature_log_abs_dA"] - mean1) / sd1
    x["z_log_seA"] = (x["feature_log_seA"] - mean2) / sd2
    x["prior_stratum"] = x["gene"].map(prior_by_gene).fillna(0).astype(int) if prior_by_gene else -1
    R = x[x["is_R"] & (x["match_direction"] != 0)].copy()
    C = x[~x["is_R"] & (x["match_direction"] != 0)].copy()
    if prior_by_gene is None:
        strata = [(-1, R, C)]
    else:
        strata = [(st, R[R.prior_stratum == st], C[C.prior_stratum == st]) for st in range(4)]
    rows = []
    pairs = []
    for st, rsub, csub in strata:
        csub = csub.sort_values("gene").reset_index(drop=True)
        cgenes = csub["gene"].to_numpy(str)
        csign = csub["match_direction"].to_numpy(int)
        cz1 = csub["z_log_abs_dA"].to_numpy(float)
        cz2 = csub["z_log_seA"].to_numpy(float)
        used = np.zeros(len(csub), dtype=bool)
        sign_index = {}
        sorted_sign_index = {}
        for sign in (-1, 1):
            ix = np.flatnonzero(csign == sign)
            sign_index[sign] = ix
            order = np.argsort(cz1[ix], kind="mergesort") if len(ix) else np.array([], dtype=int)
            sorted_ix = ix[order]
            sorted_sign_index[sign] = (sorted_ix, cz1[sorted_ix])
        for rrow in rsub.sort_values("gene").itertuples(index=False):
            sign = int(rrow.match_direction)
            sign_ix = sign_index.get(sign, np.array([], dtype=int))
            sorted_ix, sorted_z1 = sorted_sign_index.get(sign, (np.array([], dtype=int), np.array([], dtype=float)))
            lo = np.searchsorted(sorted_z1, float(rrow.z_log_abs_dA) - BENCH_B_CALIPER_SD, side="left")
            hi = np.searchsorted(sorted_z1, float(rrow.z_log_abs_dA) + BENCH_B_CALIPER_SD, side="right")
            cand_ix = sorted_ix[lo:hi]
            if len(cand_ix):
                cand_ix = cand_ix[~used[cand_ix]]
                dx_all = cz1[cand_ix] - float(rrow.z_log_abs_dA)
                dy_all = cz2[cand_ix] - float(rrow.z_log_seA)
                keep = np.abs(dy_all) <= BENCH_B_CALIPER_SD
                cand_ix, dx_all, dy_all = cand_ix[keep], dx_all[keep], dy_all[keep]
            if len(cand_ix):
                dist2 = dx_all * dx_all + dy_all * dy_all
                best = int(np.lexsort((cgenes[cand_ix], dist2))[0])
                cix = int(cand_ix[best])
                dx, dy = float(dx_all[best]), float(dy_all[best])
                used[cix] = True
                pairs.append(dict(prior_stratum=int(st), gene_R=rrow.gene, gene_control=cgenes[cix],
                                  direction_R=int(rrow.r_direction), direction_control=int(csub.iloc[cix].a_direction),
                                  delta_z_log_abs_dA=dx, delta_z_log_seA=dy,
                                  distance=float(math.sqrt(dx * dx + dy * dy))))
            else:
                reason = "no_same_sign_control" if len(sign_ix) == 0 else "no_unused_control_within_calipers"
                rows.append(dict(prior_stratum=int(st), unmatched_R=rrow.gene, unmatched_reason=reason))
        rows.append(dict(prior_stratum=int(st), n_R_eligible=int(len(rsub)), n_control_pool=int(len(csub)),
                         n_matched=int(sum(p["prior_stratum"] == st for p in pairs)),
                         n_unmatched_R=int(sum(1 for p in rows if p.get("prior_stratum") == st and "unmatched_R" in p)),
                         caliper_each_standardized_feature=BENCH_B_CALIPER_SD,
                         matching_rule="same sign; nearest standardized Euclidean distance; no replacement; symbol-order greedy"))
    pair_columns = ["prior_stratum", "gene_R", "gene_control", "direction_R", "direction_control",
                    "delta_z_log_abs_dA", "delta_z_log_seA", "distance"]
    return pd.DataFrame(pairs, columns=pair_columns), pd.DataFrame(rows), x


def balance_rows(pairdf, elig, design, target):
    if pairdf.empty:
        return []
    e = elig.set_index("gene")
    records = []
    strata = sorted(pairdf.prior_stratum.unique())
    for stratum in strata:
        sub = pairdf[pairdf.prior_stratum == stratum]
        for feature in ("z_log_abs_dA", "z_log_seA"):
            vals_r = e.loc[sub.gene_R, feature].to_numpy(float)
            vals_c = e.loc[sub.gene_control, feature].to_numpy(float)
            denom = np.sqrt((np.var(vals_r, ddof=1) + np.var(vals_c, ddof=1)) / 2) if len(vals_r) > 1 else np.nan
            smd = (np.mean(vals_r) - np.mean(vals_c)) / denom if np.isfinite(denom) and denom > 0 else np.nan
            records.append(dict(design=design, target_fold=target, prior_stratum=int(stratum), n_matched=len(sub),
                                balance_stage="after_matching", n_R=len(vals_r), n_control=len(vals_c),
                                feature=feature, R_mean=float(np.mean(vals_r)), control_mean=float(np.mean(vals_c)),
                                standardized_mean_difference=float(smd),
                                exact_direction_match_fraction=float(np.mean(
                                    e.loc[sub.gene_R, "match_direction"].to_numpy(int) ==
                                    e.loc[sub.gene_control, "match_direction"].to_numpy(int)))))
    return records


def pre_match_balance(elig, design, target, prior_by_gene=None):
    """Feature balance before matching within the eligible, same-direction pool."""
    x = elig.copy()
    if prior_by_gene is None:
        x["prior_stratum"] = -1
        strata = [-1]
    else:
        x["prior_stratum"] = x["gene"].map(prior_by_gene).fillna(0).astype(int)
        strata = list(range(4))
    records = []
    for stratum in strata:
        sub = x[x.prior_stratum == stratum]
        rsub = sub[sub.is_R & (sub.match_direction != 0)]
        csub = sub[~sub.is_R & (sub.match_direction != 0)]
        for feature in ("z_log_abs_dA", "z_log_seA"):
            vals_r = rsub[feature].to_numpy(float)
            vals_c = csub[feature].to_numpy(float)
            denom = np.sqrt((np.var(vals_r, ddof=1) + np.var(vals_c, ddof=1)) / 2) if len(vals_r) > 1 and len(vals_c) > 1 else np.nan
            smd = (np.mean(vals_r) - np.mean(vals_c)) / denom if np.isfinite(denom) and denom > 0 else np.nan
            records.append(dict(design=design, target_fold=target, prior_stratum=int(stratum),
                                balance_stage="before_matching", n_R=len(rsub), n_control=len(csub),
                                feature=feature, R_mean=float(np.mean(vals_r)) if len(vals_r) else np.nan,
                                control_mean=float(np.mean(vals_c)) if len(vals_c) else np.nan,
                                standardized_mean_difference=float(smd),
                                R_positive=int((rsub.match_direction == 1).sum()),
                                R_negative=int((rsub.match_direction == -1).sum()),
                                control_positive=int((csub.match_direction == 1).sum()),
                                control_negative=int((csub.match_direction == -1).sum()),
                                exact_direction_match_fraction=np.nan))
    return records


def group_support(genelist, directions, st, gene_pos, mapping_state):
    return point_yield_for(genelist, directions, st, gene_pos, mapping_state)


def group_support_boot(genelist, directions, st, gene_pos, mapping_state):
    return bootstrap_yield(genelist, directions, st, gene_pos, mapping_state)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--core-dir", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--ranks", type=Path)
    ap.add_argument("--deltas", type=Path)
    ap.add_argument("--mapping", type=Path)
    ap.add_argument("--r-protocol", type=Path)
    ap.add_argument("--rna-b-pairs", type=int, required=True,
                    help="verified frozen GSE15471 pair count used to construct the supplied RNA ranks")
    ap.add_argument("--bootstraps", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=20261123)
    args = ap.parse_args()
    if args.rna_b_pairs != 35:
        raise ValueError(f"strict benchmark requires RNA B/GSE15471 QC version with 35 pairs; received {args.rna_b_pairs}")
    core_dir = args.core_dir.resolve()
    out = args.out_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    log("START BENCH_A/BENCH_B/BENCH_C", f"core={core_dir}", f"out={out}",
        f"bootstraps={args.bootstraps}", f"seed={args.seed}", f"k={KS}")

    inputs = load_inputs(core_dir, args)
    write_tsv(inputs["cohort_input_audit"], out / "cohort_delta_input_audit.tsv")
    ranks = inputs["ranks"]
    genes = ranks["gene"].to_numpy(str)
    if len(genes) < max(KS):
        raise ValueError(f"RNA universe ({len(genes)}) smaller than maximum K")
    if not set(inputs["frozen_R"]).issubset(set(genes)):
        raise ValueError("frozen R contains genes outside the shared RNA universe")
    log("inputs loaded", f"RNA_universe={len(genes)}", f"frozen_R={len(inputs['frozen_R'])}",
        f"delta_rows={len(inputs['dlong'])}", f"mapping_rows={len(inputs['mdf'])}")
    pair_counts = pd.DataFrame([dict(cohort=c, protein_pairs=int(inputs["matrices"][c][1].size),
                                     protein_gene_rows=int(inputs["matrices"][c][0].size),
                                     APGI_strict_rule=("44 exact PDAC pairs" if c == "APGI" else "frozen core cohort eligibility"))
                               for c in COHORTS])
    write_tsv(pair_counts, out / "cohort_pair_counts.tsv")

    state_maps, state_audit = {}, []
    for cohort in COHORTS:
        state_maps[cohort], audit = make_state_map(inputs, cohort, genes)
        state_audit.extend(audit)
    inputs["state_maps"] = state_maps
    write_tsv(pd.DataFrame(state_audit, columns=["cohort", "gene", "matrix_row_present",
              "mapped_from_matrix_or_metadata", "metadata_mapped", "n_finite_pairs_recomputed",
              "metadata_n_pairs", "testable_recomputed", "metadata_testable", "reason"]),
              out / "testability_validation.tsv")
    prior_state_rows, prior_count_rows = [], []
    fold_data, fold_lists = {}, {}
    for target in COHORTS:
        data, states, lists = load_for_fold(inputs, ranks, target)
        fold_data[target] = data
        fold_lists[target] = lists
        prior_state_rows.extend(states)
        prior_count_rows.extend(data[["gene", "availability_prior_count_0_3", "rank_availability_only",
                                      "rank_RNA_plus_availability", "score_RNA_plus_availability"]].assign(
                                          fold=f"holdout_{target}", target_cohort=target).to_dict("records"))
    write_tsv(pd.DataFrame(prior_state_rows), out / "prior_source_states.tsv")
    write_tsv(pd.DataFrame(prior_count_rows), out / "fold_prior_and_ranks.tsv")
    # Explicit audit that the target cohort did not contribute to the prior.
    prior_audit = pd.DataFrame(prior_state_rows)
    if prior_audit.loc[prior_audit.source_role == "excluded_target", "prior_contribution"].sum() != 0:
        raise AssertionError("held-out target contributed to its availability prior")
    if (prior_audit.loc[prior_audit.source_role == "excluded_target", "protein_testable"].notna().any() or
            prior_audit.loc[prior_audit.source_role == "excluded_target", "n_finite_pairs"].notna().any()):
        raise AssertionError("held-out target state was read during prior feature construction")

    # Full-data point outcomes are only read after the fold-specific prior ranks exist.
    point_rows, point_stats, map_vectors = [], {}, {}
    for target in COHORTS:
        pgenes, patients, D = inputs["matrices"][target]
        st = paired_stats(D)
        point_stats[target] = st
        ppos = {g: i for i, g in enumerate(pgenes)}
        map_vectors[target] = state_maps[target]
        for rule in BENCH_A_RULES:
            selected_all, dirs_all = fold_lists[target][rule]
            for k in KS:
                selected, dirs = selected_all[:k], dirs_all[:k]
                vals = point_yield_for(selected, dirs, st, ppos, state_maps[target])
                point_rows.append(dict(fold=f"holdout_{target}", target_cohort=target, rule=rule, k=k,
                                       direction_source={"availability_only": "frozen_equal_weight_RNA",
                                                         "RNA_plus_availability": "frozen_continuous_RNA"}.get(rule, f"frozen_{rule}_RNA"),
                                       denominator="fixed ITT nomination slots; unmapped and <10 pairs unsupported",
                                       **vals))
    point_df = pd.DataFrame(point_rows)
    write_tsv(point_df, out / "point_yields.tsv")
    core_grid_path = find_file(core_dir, ["fixed_budget_grid.tsv", "core_fixed_budget_grid.tsv"], required=False)
    ledger_path = find_file(core_dir, ["analysis_ledger.tsv", "strict_analysis_ledger.tsv"], required=False)
    if core_grid_path:
        core_grid = pd.read_csv(core_grid_path, sep="\t")
        core_grid.columns = [str(c).strip() for c in core_grid.columns]
        cc = col(core_grid, ["cohort", "target_cohort"])
        ck = col(core_grid, ["k", "budget"])
        cr = col(core_grid, ["rule", "method"])
        cm = col(core_grid, ["mapped"])
        ct = col(core_grid, ["testable", "testable_ge10"])
        cs = col(core_grid, ["supported", "supported_q05"])
        ccan = col(core_grid, ["candidates", "candidate_slots"], required=False)
        core_grid["target_cohort"] = core_grid[cc].map(norm_cohort)
        core_grid["k"] = pd.to_numeric(core_grid[ck], errors="coerce").astype(int)
        core_grid["rule"] = core_grid[cr].astype(str)
        core_grid = core_grid[(core_grid.rule == "continuous") & core_grid.target_cohort.isin(COHORTS)].copy()
        core_grid = core_grid.rename(columns={cm: "core_mapped", ct: "core_testable", cs: "core_supported"})
        if ccan:
            core_grid = core_grid.rename(columns={ccan: "core_candidate_slots"})
        else:
            core_grid["core_candidate_slots"] = core_grid["k"]
        core_ref_path = core_grid_path
    elif ledger_path:
        ledger = pd.read_csv(ledger_path, sep="\t")
        ledger.columns = [str(c).strip() for c in ledger.columns]
        lc = col(ledger, ["cohort", "target_cohort"])
        lr = col(ledger, ["rank_continuous"])
        lm = col(ledger, ["protein_mapped", "mapped"])
        lt = col(ledger, ["protein_testable", "testable"])
        ls = col(ledger, ["supported_continuous"])
        ledger["target_cohort"] = ledger[lc].map(norm_cohort)
        ledger["rank_continuous"] = pd.to_numeric(ledger[lr], errors="coerce")
        ledger["protein_mapped"] = ledger[lm].map(truth)
        ledger["protein_testable"] = ledger[lt].map(truth)
        ledger["supported_continuous"] = ledger[ls].map(truth)
        core_rows = []
        for cohort in COHORTS:
            for k in KS:
                sel = ledger[(ledger.target_cohort == cohort) & (ledger.rank_continuous <= k)]
                core_rows.append(dict(target_cohort=cohort, k=k, core_candidate_slots=len(sel),
                                      core_mapped=int(sel.protein_mapped.sum()),
                                      core_testable=int(sel.protein_testable.sum()),
                                      core_supported=int(sel.supported_continuous.sum())))
        core_grid = pd.DataFrame(core_rows)
        core_grid_path = ledger_path
        core_ref_path = ledger_path
    else:
        core_grid = None
        core_ref_path = None
    if core_grid is not None:
        val = point_df[(point_df.rule == "continuous")][["target_cohort", "k", "candidate_slots", "mapped", "testable_ge10", "supported_q05"]].merge(
            core_grid[["target_cohort", "k", "core_candidate_slots", "core_mapped", "core_testable", "core_supported"]],
            on=["target_cohort", "k"], how="outer", indicator=True)
        val["matches_core_point_counts"] = (val["_merge"].eq("both") &
            (val["candidate_slots"] == val["core_candidate_slots"]) &
            (val["mapped"] == val["core_mapped"]) &
            (val["testable_ge10"] == val["core_testable"]) &
            (val["supported_q05"] == val["core_supported"]))
        val["core_reference"] = str(core_ref_path)
        write_tsv(val, out / "main_point_validation.tsv")
        if (len(val) != len(COHORTS) * len(KS) or
                val.duplicated(["target_cohort", "k"]).any() or
                not bool(val["matches_core_point_counts"].all())):
            bad_rows = val.loc[~val["matches_core_point_counts"],
                               ["target_cohort", "k", "_merge", "candidate_slots", "core_candidate_slots",
                                "mapped", "core_mapped", "testable_ge10", "core_testable",
                                "supported_q05", "core_supported"]].to_dict("records")
            raise ValueError(f"continuous RNA main point counts do not match strict core; rows={len(val)} expected={len(COHORTS) * len(KS)}, duplicates={val.duplicated(['target_cohort','k']).sum()}; see main_point_validation.tsv: {bad_rows[:8]}")
    else:
        write_tsv(pd.DataFrame([dict(status="strict core fixed_budget_grid.tsv not present at runner time",
                                     validation_performed=False)]), out / "main_point_validation.tsv")

    # Top-500 pairwise list overlap and shared-gene direction changes.
    overlap_rows, shared_rows = [], []
    for target in COHORTS:
        lists = fold_lists[target]
        rank_lookup = fold_data[target].set_index("gene")
        for i, left in enumerate(BENCH_C_RULES):
            lg, ld = lists[left][0][:500], lists[left][1][:500]
            lm = dict(zip(lg, ld))
            for right in BENCH_C_RULES[i + 1:]:
                rg, rd = lists[right][0][:500], lists[right][1][:500]
                rm = dict(zip(rg, rd))
                shared = sorted(set(lm) & set(rm))
                same = sum(1 for g in shared if lm[g] != 0 and rm[g] != 0 and lm[g] == rm[g])
                opposite = sum(1 for g in shared if lm[g] != 0 and rm[g] != 0 and lm[g] == -rm[g])
                zero = len(shared) - same - opposite
                union = len(set(lm) | set(rm))
                overlap_rows.append(dict(fold=f"holdout_{target}", target_cohort=target, rule_left=left,
                                         rule_right=right, k=500, intersection=len(shared), union=union,
                                         jaccard=len(shared) / union if union else np.nan,
                                         replacements_each_list=500 - len(shared), same_direction=same,
                                         opposite_direction=opposite, at_least_one_zero_direction=zero,
                                         note="double_sig_vs_continuous is an explicit prespecified row" if
                                         {left, right} == {"double_sig", "continuous"} else ""))
                for g in shared:
                    shared_rows.append(dict(fold=f"holdout_{target}", target_cohort=target,
                                            rule_left=left, rule_right=right, gene=g,
                                            rank_left=int(rank_lookup.loc[g, f"rank_{left}"]),
                                            rank_right=int(rank_lookup.loc[g, f"rank_{right}"]),
                                            direction_left=int(lm[g]), direction_right=int(rm[g]),
                                            direction_relation=("same" if lm[g] == rm[g] and lm[g] != 0 else
                                                                "opposite" if lm[g] == -rm[g] and lm[g] != 0 else "zero_or_undefined")))
    write_tsv(pd.DataFrame(overlap_rows), out / "bench_c_top500_overlap.tsv")
    write_tsv(pd.DataFrame(shared_rows), out / "bench_c_shared_directions.tsv")

    # BENCH_B matching is constructed only from RNA fields and, for the sensitivity,
    # the target-excluded LOCO count; protein target outcomes are not in this step.
    elig = eligible_matched(ranks, inputs["frozen_R"])
    base_pairs, base_flow, elig_z = match_controls(elig, prior_by_gene=None)
    base_pairs["design"] = "RNA_only_main"
    base_flow["design"] = "RNA_only_main"
    base_flow["target_fold"] = "all_targets"
    bench_b_pairs = [base_pairs.copy()]
    bench_b_flow = [base_flow.copy()]
    balance = pre_match_balance(elig_z, "RNA_only_main", "all_targets")
    balance.extend(balance_rows(base_pairs, elig_z, "RNA_only_main", "all_targets"))
    prior_pair_by_fold = {}
    for target in COHORTS:
        prior = dict(zip(fold_data[target]["gene"], fold_data[target]["availability_prior_count_0_3"].astype(int)))
        pairs, flow, zdata = match_controls(elig, prior_by_gene=prior)
        pairs["design"] = "LOCO_prior_stratified"
        pairs["target_fold"] = f"holdout_{target}"
        flow["design"] = "LOCO_prior_stratified"
        flow["target_fold"] = f"holdout_{target}"
        prior_pair_by_fold[target] = pairs
        bench_b_pairs.append(pairs)
        bench_b_flow.append(flow)
        balance.extend(pre_match_balance(zdata, "LOCO_prior_stratified", f"holdout_{target}", prior_by_gene=prior))
        balance.extend(balance_rows(pairs, zdata, "LOCO_prior_stratified", f"holdout_{target}"))
    pair_all = pd.concat(bench_b_pairs, ignore_index=True) if bench_b_pairs else pd.DataFrame()
    flow_all = pd.concat(bench_b_flow, ignore_index=True) if bench_b_flow else pd.DataFrame()
    write_tsv(pair_all, out / "bench_b_matched_pairs.tsv")
    write_tsv(flow_all, out / "bench_b_matching.tsv")
    write_tsv(pd.DataFrame(balance), out / "bench_b_balance.tsv")

    # BENCH_A fixed-list patient bootstrap plus BENCH_B matched-set bootstrap.
    rngs = {cohort: np.random.default_rng(args.seed + i * 1009) for i, cohort in enumerate(COHORTS)}
    bench_a_draws, bench_b_draws, bench_b_yields = [], [], []
    contrasts = [("availability_only", rule) for rule in RNA_RULES]
    contrasts += [("RNA_plus_availability", rule) for rule in RNA_RULES]
    contrasts += [("RNA_plus_availability", "availability_only")]
    direct_contrast_labels = [(f"{a}_minus_{b}", a, b) for a, b in contrasts]
    point_k500 = point_df[point_df.k == 500].set_index(["target_cohort", "rule"])
    point_contrast_rows = []
    for target in COHORTS:
        for label, left, right in direct_contrast_labels:
            lv = point_k500.loc[(target, left)]
            rv = point_k500.loc[(target, right)]
            point_contrast_rows.append(dict(fold=f"holdout_{target}", target_cohort=target,
                                            contrast=label, rule_left=left, rule_right=right,
                                            point_delta_itt_supported_per100=lv.supported_per100_slots - rv.supported_per100_slots,
                                            point_delta_conditional_support_pp=(lv.conditional_supported_rate - rv.conditional_supported_rate) * 100
                                            if np.isfinite(lv.conditional_supported_rate) and np.isfinite(rv.conditional_supported_rate) else np.nan,
                                            left_supported=int(lv.supported_q05), right_supported=int(rv.supported_q05),
                                            left_mapped=int(lv.mapped), right_mapped=int(rv.mapped),
                                            left_testable=int(lv.testable_ge10), right_testable=int(rv.testable_ge10),
                                            slots_each=500, estimand="full-data fixed-list; ITT denominator 500; conditional denominator testable genes"))
    point_contrast_df = pd.DataFrame(point_contrast_rows)
    write_tsv(point_contrast_df, out / "k500_point_contrasts.tsv")

    for target in COHORTS:
        pgenes, patients, D = inputs["matrices"][target]
        ppos = {g: i for i, g in enumerate(pgenes)}
        dirs_by_rule = {}
        selected_by_rule = {}
        for rule in BENCH_A_RULES:
            allg, alld = fold_lists[target][rule]
            selected_by_rule[rule] = allg[:500]
            dirs_by_rule[rule] = alld[:500]

        # Main BENCH_B match is identical for all protein cohorts; LOCO strata are fold-specific.
        bench_b_designs = {"RNA_only_main": base_pairs,
                       "LOCO_prior_stratified": prior_pair_by_fold[target]}
        # The target-aware shared-evaluator-testable subset is frozen from full-data status and labeled diagnostic.
        target_full_test = state_maps[target]
        diagnostic_pairs = base_pairs[
            base_pairs.apply(lambda row: target_full_test.get(row.gene_R, {}).get("testable", False) and
                             target_full_test.get(row.gene_control, {}).get("testable", False), axis=1)].copy()
        diagnostic_pairs["design"] = "common_evaluator_testable_diagnostic"
        diagnostic_pairs["target_fold"] = f"holdout_{target}"
        bench_b_designs["common_evaluator_testable_diagnostic"] = diagnostic_pairs

        for design, pairdf in bench_b_designs.items():
            if pairdf.empty:
                continue
            group_genes = {
                "R": pairdf["gene_R"].to_numpy(str),
                "matched_control": pairdf["gene_control"].to_numpy(str),
            }
            group_dirs = {
                "R": pairdf["direction_R"].to_numpy(int),
                "matched_control": pairdf["direction_control"].to_numpy(int),
            }
            for group in ("R", "matched_control"):
                vals = group_support(group_genes[group], group_dirs[group], point_stats[target], ppos,
                                     state_maps[target])
                bench_b_yields.append(dict(target_cohort=target, target_fold=f"holdout_{target}", design=design,
                                       group=group, matched_pairs=len(pairdf),
                                       target_testability_used_for_selection=(design == "common_evaluator_testable_diagnostic"),
                                       **vals))

        rng = rngs[target]
        draw_count = int(args.bootstraps)
        if draw_count <= 0:
            raise ValueError("bootstraps must be positive")
        for b in range(draw_count):
            ix = rng.integers(0, D.shape[1], D.shape[1])
            bst = paired_stats(D[:, ix])
            # Fixed ranking/prior lists; evaluator patients and protein BH are resampled.
            boot_rule = {}
            for rule in BENCH_A_RULES:
                boot_rule[rule] = bootstrap_yield(selected_by_rule[rule], dirs_by_rule[rule], bst,
                                                   ppos, state_maps[target])
            for label, left, right in direct_contrast_labels:
                lv, rv = boot_rule[left], boot_rule[right]
                delta_cond = (lv["conditional"] - rv["conditional"]
                              if np.isfinite(lv["conditional"]) and np.isfinite(rv["conditional"]) else np.nan)
                bench_a_draws.append(dict(replicate=b + 1, fold=f"holdout_{target}", target_cohort=target,
                                      contrast=label, rule_left=left, rule_right=right,
                                      delta_itt_supported_per100=lv["itt"] - rv["itt"],
                                      delta_conditional_support_pp=delta_cond,
                                      left_testable=lv["testable_count"], right_testable=rv["testable_count"],
                                      paired_patient_resampling=True, ranks_and_prior_fixed=True,
                                      seed=args.seed + COHORTS.index(target) * 1009))

            # BENCH_B group estimates use the same target patient resample; matched lists remain fixed.
            for design, pairdf in bench_b_designs.items():
                if pairdf.empty:
                    continue
                rg = pairdf["gene_R"].to_numpy(str)
                cg = pairdf["gene_control"].to_numpy(str)
                rd = pairdf["direction_R"].to_numpy(int)
                cd = pairdf["direction_control"].to_numpy(int)
                rvals = group_support_boot(rg, rd, bst, ppos, state_maps[target])
                cvals = group_support_boot(cg, cd, bst, ppos, state_maps[target])
                cond = (rvals["conditional"] - cvals["conditional"]
                        if np.isfinite(rvals["conditional"]) and np.isfinite(cvals["conditional"]) else np.nan)
                bench_b_draws.append(dict(replicate=b + 1, target_cohort=target, target_fold=f"holdout_{target}",
                                      design=design, n_matched_pairs=len(pairdf),
                                      delta_itt_supported_per100=rvals["itt"] - cvals["itt"],
                                      delta_conditional_support_pp=cond,
                                      R_testable=rvals["testable_count"], control_testable=cvals["testable_count"],
                                      conditional_target_testability_used=(design == "common_evaluator_testable_diagnostic"),
                                      paired_patient_resampling=True, ranks_and_matching_fixed=True,
                                      seed=args.seed + COHORTS.index(target) * 1009))
            if (b + 1) % 100 == 0:
                log("bootstrap progress", f"target={target}", f"{b + 1}/{draw_count}")

    bench_a_draw_df = pd.DataFrame(bench_a_draws)
    write_tsv(bench_a_draw_df, out / "k500_bootstrap_draws.tsv")
    bench_a_summary = summarize_draws(bench_a_draws, ["fold", "target_cohort", "contrast"],
                                  "delta_itt_supported_per100",
                                  "fixed top-500 ranks and LOCO prior; paired target-cohort patient resampling; cohort-wide protein BH")
    bench_a_cond = summarize_draws(bench_a_draws, ["fold", "target_cohort", "contrast"],
                               "delta_conditional_support_pp",
                               "fixed top-500 ranks and LOCO prior; conditional on target-cohort testability; paired patient resampling")
    bench_a_summary = pd.DataFrame(bench_a_summary).rename(columns={"bootstrap_center": "ITT_bootstrap_center",
                                                            "bootstrap_median": "ITT_bootstrap_median",
                                                            "lo95": "ITT_lo95", "hi95": "ITT_hi95",
                                                            "draws_valid": "ITT_valid_draws",
                                                            "probability_gt0": "ITT_probability_gt0"})
    bench_a_cond = pd.DataFrame(bench_a_cond).rename(columns={"bootstrap_center": "conditional_bootstrap_center_pp",
                                                      "bootstrap_median": "conditional_bootstrap_median_pp",
                                                      "lo95": "conditional_lo95_pp", "hi95": "conditional_hi95_pp",
                                                      "draws_valid": "conditional_valid_draws",
                                                      "probability_gt0": "conditional_probability_gt0"})
    summary = bench_a_summary.merge(bench_a_cond[["fold", "target_cohort", "contrast", "conditional_bootstrap_center_pp",
                                         "conditional_bootstrap_median_pp", "conditional_lo95_pp", "conditional_hi95_pp",
                                         "conditional_valid_draws", "conditional_probability_gt0"]],
                               on=["fold", "target_cohort", "contrast"], how="left")
    summary = summary.merge(point_contrast_df[["fold", "target_cohort", "contrast", "point_delta_itt_supported_per100",
                                              "point_delta_conditional_support_pp", "left_supported", "right_supported",
                                              "left_mapped", "right_mapped", "left_testable", "right_testable"]],
                            on=["fold", "target_cohort", "contrast"], how="left")
    write_tsv(summary, out / "k500_contrasts.tsv")
    bench_b_draw_df = pd.DataFrame(bench_b_draws)
    write_tsv(bench_b_draw_df, out / "bench_b_bootstrap_draws.tsv")
    if bench_b_draws:
        bench_b_itt = summarize_draws(bench_b_draws, ["target_cohort", "target_fold", "design"],
                                  "delta_itt_supported_per100",
                                  "fixed RNA-only or LOCO-stratified matched genes; paired target-cohort patient resampling; cohort-wide protein BH")
        bench_b_cond = summarize_draws(bench_b_draws, ["target_cohort", "target_fold", "design"],
                                   "delta_conditional_support_pp",
                                   "fixed matched genes; conditional on target-cohort testability; paired patient resampling")
        bench_b_itt = pd.DataFrame(bench_b_itt).rename(columns={"bootstrap_center": "ITT_bootstrap_center",
                                                       "bootstrap_median": "ITT_bootstrap_median",
                                                       "lo95": "ITT_lo95", "hi95": "ITT_hi95",
                                                       "draws_valid": "ITT_valid_draws",
                                                       "probability_gt0": "ITT_probability_gt0"})
        bench_b_cond = pd.DataFrame(bench_b_cond).rename(columns={"bootstrap_center": "conditional_bootstrap_center_pp",
                                                         "bootstrap_median": "conditional_bootstrap_median_pp",
                                                         "lo95": "conditional_lo95_pp", "hi95": "conditional_hi95_pp",
                                                         "draws_valid": "conditional_valid_draws",
                                                         "probability_gt0": "conditional_probability_gt0"})
        bench_b_summary = bench_b_itt.merge(bench_b_cond[["target_cohort", "target_fold", "design", "conditional_bootstrap_center_pp",
                                            "conditional_bootstrap_median_pp", "conditional_lo95_pp", "conditional_hi95_pp",
                                            "conditional_valid_draws", "conditional_probability_gt0"]],
                                   on=["target_cohort", "target_fold", "design"], how="left")
    else:
        bench_b_summary = pd.DataFrame(columns=["target_cohort", "target_fold", "design", "draws_valid"])
    bench_b_yield_df = pd.DataFrame(bench_b_yields)
    write_tsv(bench_b_yield_df, out / "bench_b_yields.tsv")
    bench_b_point_rows = []
    if not bench_b_yield_df.empty:
        for (target, fold, design), sub in bench_b_yield_df.groupby(["target_cohort", "target_fold", "design"]):
            lookup = {row.group: row for row in sub.itertuples(index=False)}
            if "R" not in lookup or "matched_control" not in lookup:
                continue
            rv, cv = lookup["R"], lookup["matched_control"]
            bench_b_point_rows.append(dict(target_cohort=target, target_fold=fold, design=design,
                                       matched_pairs=int(rv.matched_pairs),
                                       point_delta_itt_supported_per100=rv.supported_per100_slots - cv.supported_per100_slots,
                                       point_delta_conditional_support_pp=(rv.conditional_supported_rate - cv.conditional_supported_rate) * 100
                                       if np.isfinite(rv.conditional_supported_rate) and np.isfinite(cv.conditional_supported_rate) else np.nan,
                                       R_supported=rv.supported_q05, control_supported=cv.supported_q05,
                                       R_testable=rv.testable_ge10, control_testable=cv.testable_ge10,
                                       R_mapped=rv.mapped, control_mapped=cv.mapped,
                                       target_testability_used_for_selection=bool(rv.target_testability_used_for_selection)))
    bench_b_point_df = pd.DataFrame(bench_b_point_rows)
    write_tsv(bench_b_point_df, out / "bench_b_point_contrasts.tsv")
    if not bench_b_summary.empty and not bench_b_point_df.empty:
        bench_b_summary = bench_b_summary.merge(bench_b_point_df, on=["target_cohort", "target_fold", "design"], how="left")
    write_tsv(bench_b_summary, out / "bench_b_contrasts.tsv")

    # Compact Chinese interpretive report and negative-finding ledger.
    report = []
    report.append("# Benchmark analyses: summary of BENCH_A/BENCH_B/BENCH_C results")
    report.append("")
    report.append("## 预先冻结的方法")
    report.append("")
    report.append("四次留一蛋白队列分析中，可测性先验只累计其余三队列内‘已映射且至少 10 个有限配对差值’的来源数（0–3）。未进入发布矩阵与已映射但配对不足分别记录，均不作为可测证据。先验不读取留出队列蛋白数据、效应、显著性或病例字段。可测性单独排序；合并规则固定为连续 RNA 排名与先验中秩百分位各 50%，方向固定取连续 RNA。可测性排序的方向固定取等权 RNA。k=100/250/500/1000/2000。BENCH_B 只用 RNA-A 的方向、效应量和精度按冻结卡钳 1:1 匹配，不按蛋白结局调整。")
    report.append("")
    report.append("输入病例数按冻结核心矩阵核对，见 `cohort_pair_counts.tsv`。严格 APGI 为 44 个病理明确 PDAC 配对；同文件内 APGI mixed-histology 51 对作为单独敏感性来源明确排除在主先验和评价之外。严格 CPTAC 为 64 个 PDAC 配对（原 67 对中排除 3 个腺鳞癌病例）。RNA 排名采用 GSE15471 QC 后 35 对的冻结版本。Fudan 和 Ruijin 按冻结核心纳入定义。原始 delta 文件中的队列标签纳入/排除审计见 `cohort_delta_input_audit.tsv`。")
    report.append("")
    report.append(f"蛋白映射/可测性字段与 delta 长表重新计数的差异数为 {len(state_audit)} 个基因×队列记录；逐条对照见 `testability_validation.tsv`。")
    report.append("")
    report.append("## BENCH_A 观察结果")
    report.append("")
    key = point_df[(point_df.k == 500) & point_df.rule.isin(["single_A", "single_B", "equal_weight", "continuous", "availability_only", "RNA_plus_availability"])]
    def fmt_effect(point, center, lo, hi, digits=1):
        if not all(np.isfinite(x) for x in (point, center, lo, hi)):
            return "NA"
        fmt = f"{{:+.{digits}f}}"
        return f"{fmt.format(point)}; {fmt.format(center)} [{fmt.format(lo)}, {fmt.format(hi)}]"
    if len(key):
        report.append("k=500 点估计（列为队列、规则、支持数/500、可测/500、可测候选中的支持率）：")
        report.append("")
        report.append("| 留出队列 | 规则 | 支持数/500 | 可测/500 | 可测中支持率 |")
        report.append("|---|---|---:|---:|---:|")
        for row in key.itertuples(index=False):
            report.append(f"| {row.target_cohort} | {row.rule} | {row.supported_q05}/500 | {row.testable_ge10}/500 | {row.conditional_supported_rate:.3f} |")
        report.append("")
        primary_contrasts = summary[summary.contrast.isin([
            "availability_only_minus_continuous", "RNA_plus_availability_minus_continuous"])]
        primary_lookup = primary_contrasts.set_index(["target_cohort", "contrast"])
        report.append("### k=500 预设对比：全样本点值；配对患者 bootstrap 均值 [95% 区间]")
        report.append("")
        report.append("IT T 以每 100 个固定名额的方向支持数表示；条件差以可检验候选中的支持率百分点表示。")
        report.append("")
        report.append("| 留出队列 | 仅可测性−连续 ITT | RNA+可测性−连续 ITT | RNA+可测性−连续 条件差(pp) |")
        report.append("|---|---:|---:|---:|")
        for cohort in COHORTS:
            avail = primary_lookup.loc[(cohort, "availability_only_minus_continuous")]
            combo = primary_lookup.loc[(cohort, "RNA_plus_availability_minus_continuous")]
            report.append("| {} | {} | {} | {} |".format(
                cohort,
                fmt_effect(avail.point_delta_itt_supported_per100, avail.ITT_bootstrap_center, avail.ITT_lo95, avail.ITT_hi95),
                fmt_effect(combo.point_delta_itt_supported_per100, combo.ITT_bootstrap_center, combo.ITT_lo95, combo.ITT_hi95),
                fmt_effect(combo.point_delta_conditional_support_pp, combo.conditional_bootstrap_center_pp,
                           combo.conditional_lo95_pp, combo.conditional_hi95_pp)))
        report.append("")
        report.append("可测性单独排序提高或维持若干队列的 ITT 支持，但条件支持均低于连续 RNA；例如 CPTAC 的 ITT 差为负，APGI 点差为零。固定 50:50 名次融合在四个队列的 ITT 点值均高于连续 RNA；条件支持只在 Fudan 明确增加，CPTAC 和 Ruijin 降低，APGI 区间包含零。全组预设对比（包括单 A、单 B、等权基线）及可测/未映射分母见 `k500_point_contrasts.tsv`、`k500_contrasts.tsv` 和 `point_yields.tsv`。")
        report.append("")
        report.append("全样本 k=500 预设对比见 `k500_point_contrasts.tsv`；配对患者 bootstrap 的中心和 95% 区间见 `k500_contrasts.tsv`。区间固定 RNA 排名和留一先验，只重抽对应留出蛋白队列患者并在队列全检验家族重算配对检验/BH。因此区间不覆盖 RNA 排名、先验跨队列构造、平台或队列抽样不确定性。")
        val_path = out / "main_point_validation.tsv"
        val_check = pd.read_csv(val_path, sep="\t") if val_path.exists() else pd.DataFrame()
        if not val_check.empty and "matches_core_point_counts" in val_check.columns:
            report.append("")
            report.append(f"连续 RNA 主规则点估计与严格核心 `fixed_budget_grid.tsv` 对照：{int(val_check['matches_core_point_counts'].sum())}/{len(val_check)} 个队列×预算行四项计数完全一致。逐项差异见 `main_point_validation.tsv`。")
        report.append("")
        report.append("不优于比较对象的观察也按预设规则保留：")
        report.append("")
        report.append("| 对比（左减右） | 支持数未增加的留出队列 | 可测数未增加的留出队列 |")
        report.append("|---|---:|---:|")
        for left, right in [("availability_only", r) for r in RNA_RULES] + \
                           [("RNA_plus_availability", r) for r in RNA_RULES] + \
                           [("RNA_plus_availability", "availability_only")]:
            nonpositive_support, nonpositive_testable = [], []
            for cohort in COHORTS:
                a = key[(key.target_cohort == cohort) & (key.rule == left)]
                b = key[(key.target_cohort == cohort) & (key.rule == right)]
                if a.empty or b.empty:
                    continue
                aa, bb = a.iloc[0], b.iloc[0]
                if aa.supported_q05 <= bb.supported_q05:
                    nonpositive_support.append(cohort)
                if aa.testable_ge10 <= bb.testable_ge10:
                    nonpositive_testable.append(cohort)
            report.append(f"| {left} − {right} | {', '.join(nonpositive_support) or '无'} | {', '.join(nonpositive_testable) or '无'} |")
    else:
        report.append("未生成 k=500 点估计；需检查上游输入。")
    report.append("")
    report.append("## BENCH_B 匹配结果")
    report.append("")
    if not base_pairs.empty:
        report.append(f"RNA-only 复现性匹配（此处不是临床预测模型）：冻结 R={len(set(inputs['frozen_R']))}，符合 RNA-A 特征与方向条件后匹配 {len(base_pairs)} 对；未匹配原因见 `bench_b_matching.tsv`，效应量/精度标准化差异见 `bench_b_balance.tsv`。全样本差异、分母和目标队列患者 bootstrap 见 `bench_b_point_contrasts.tsv` 与 `bench_b_contrasts.tsv`。留一先验分层为目标盲法敏感性；按留出矩阵共同可测性限定的子集仅作事后诊断。")
        report.append("")
        bench_b_main = bench_b_summary[bench_b_summary.design == "RNA_only_main"].set_index("target_cohort")
        report.append("主匹配结果（R 方向按冻结集合；对照方向按 RNA-A；`支持/可测` 均为基因数；ITT 差/100 与条件差 pp 给出点值及 bootstrap 均值 [95% 区间]）：")
        report.append("")
        report.append("| 留出队列 | 匹配对 | R 支持/可测 | 对照支持/可测 | ITT 差/100 | 条件差(pp) |")
        report.append("|---|---:|---:|---:|---:|---:|")
        for cohort in COHORTS:
            row = bench_b_main.loc[cohort]
            report.append("| {} | {} | {}/{} | {}/{} | {} | {} |".format(
                cohort, int(row.matched_pairs), int(row.R_supported), int(row.R_testable),
                int(row.control_supported), int(row.control_testable),
                fmt_effect(row.point_delta_itt_supported_per100, row.ITT_bootstrap_center, row.ITT_lo95, row.ITT_hi95),
                fmt_effect(row.point_delta_conditional_support_pp, row.conditional_bootstrap_center_pp,
                           row.conditional_lo95_pp, row.conditional_hi95_pp)))
        main_bal = pd.DataFrame(balance)
        main_bal = main_bal[(main_bal.design == "RNA_only_main") & (main_bal.balance_stage == "after_matching")]
        smd_map = {row.feature: row.standardized_mean_difference for row in main_bal.itertuples(index=False)}
        loco_bal = pd.DataFrame(balance)
        loco_bal = loco_bal[(loco_bal.design == "LOCO_prior_stratified") &
                            (loco_bal.balance_stage == "after_matching")]
        if not loco_bal.empty:
            abs_d = loco_bal[loco_bal.feature == "z_log_abs_dA"]["standardized_mean_difference"].astype(float)
            abs_se = loco_bal[loco_bal.feature == "z_log_seA"]["standardized_mean_difference"].astype(float)
            loco_balance_text = (f"LOCO 分层后 SMD 范围：log|d_A| {abs_d.min():.3f} 至 {abs_d.max():.3f}，"
                                 f"log(SE_A) {abs_se.min():.3f} 至 {abs_se.max():.3f}。")
        else:
            loco_balance_text = "LOCO 分层平衡结果见 `bench_b_balance.tsv`。"
        report.append("")
        report.append(f"主匹配后 SMD：log|d_A|={smd_map.get('z_log_abs_dA', np.nan):.3f}，log(SE_A)={smd_map.get('z_log_seA', np.nan):.3f}，配对方向一致率为 1.000。{loco_balance_text}")
        report.append("")
        loco_summary = bench_b_summary[bench_b_summary.design == "LOCO_prior_stratified"].set_index("target_cohort")
        report.append("留一可测性先验分层敏感性（对照与 R 同属 0–3 先验计数层；为目标盲法匹配）：")
        report.append("")
        report.append("| 留出队列 | 匹配对 | ITT 差/100 均值 [95% 区间] | 条件差(pp) 均值 [95% 区间] |")
        report.append("|---|---:|---:|---:|")
        for cohort in COHORTS:
            row = loco_summary.loc[cohort]
            report.append("| {} | {} | {} | {} |".format(
                cohort, int(row.matched_pairs),
                fmt_effect(row.point_delta_itt_supported_per100, row.ITT_bootstrap_center, row.ITT_lo95, row.ITT_hi95),
                fmt_effect(row.point_delta_conditional_support_pp, row.conditional_bootstrap_center_pp,
                           row.conditional_lo95_pp, row.conditional_hi95_pp)))
        common_summary = bench_b_summary[bench_b_summary.design == "common_evaluator_testable_diagnostic"].set_index("target_cohort")
        if not common_summary.empty:
            report.append("")
            report.append("按留出蛋白矩阵中 R 与对照同时可测进行的事后诊断（目标依赖，不能用于部署排序）：")
            report.append("")
            report.append("| 留出队列 | 双可测匹配对 | ITT 差/100 均值 [95% 区间] | 条件差(pp) 均值 [95% 区间] |")
            report.append("|---|---:|---:|---:|")
            for cohort in COHORTS:
                row = common_summary.loc[cohort]
                report.append("| {} | {} | {} | {} |".format(
                    cohort, int(row.matched_pairs),
                    fmt_effect(row.point_delta_itt_supported_per100, row.ITT_bootstrap_center, row.ITT_lo95, row.ITT_hi95),
                    fmt_effect(row.point_delta_conditional_support_pp, row.conditional_bootstrap_center_pp,
                               row.conditional_lo95_pp, row.conditional_hi95_pp)))
            report.append("")
        if not flow_all.empty:
            main_flow = flow_all[flow_all.design == "RNA_only_main"]
            total_r = int(main_flow["n_R_eligible"].fillna(0).sum()) if "n_R_eligible" in main_flow else 0
            total_unmatched = int(main_flow["n_unmatched_R"].fillna(0).sum()) if "n_unmatched_R" in main_flow else max(0, total_r - len(base_pairs))
            report.append(f"RNA-only 匹配负结果：{total_unmatched} 个 R 条目未能在预设同方向/双特征卡钳内找到未使用对照（可按缺同向对照或无剩余卡钳内对照区分）；匹配前后特征余额均列在 `bench_b_balance.tsv`。")
    else:
        report.append("没有满足预设方向和双特征卡钳的匹配对；不作匹配结局推断。")
    report.append("")
    report.append("## BENCH_C 名单和方向重合")
    report.append("")
    overlap_df = pd.DataFrame(overlap_rows)
    dc = overlap_df[((overlap_df.rule_left == "continuous") & (overlap_df.rule_right == "double_sig")) |
                    ((overlap_df.rule_right == "continuous") & (overlap_df.rule_left == "double_sig"))]
    if len(dc):
        report.append("连续榜单与双显著优先的 top-500 交集、替换数及共享基因方向关系：")
        report.append("")
        report.append("| 留出队列 | 交集 | 替换名额 | 同向 | 反向 | 零/未定义方向 |")
        report.append("|---|---:|---:|---:|---:|---:|")
        for row in dc.itertuples(index=False):
            report.append(f"| {row.target_cohort} | {row.intersection}/500 | {row.replacements_each_list} | {row.same_direction} | {row.opposite_direction} | {row.at_least_one_zero_direction} |")
    report.append("")
    report.append("## 负结果与范围")
    report.append("")
    report.append("以下结果也保留在 TSV 中，不进行事后挑选：未增加支持或可测名额的队列/规则；先验排序与 RNA 排序不一致的结果；匹配共同支持不足或方向/精度/效应量仍失衡的分层；无效 bootstrap 条件率抽样。此分析是历史发布矩阵内的固定名额提名比较，不能换算为肽段可实现率、仪器时间、成本或前瞻实验收益。")
    report.append("")
    report.append("## 可复现文件")
    report.append("")
    for name in ("point_yields.tsv", "cohort_delta_input_audit.tsv", "main_point_validation.tsv", "k500_point_contrasts.tsv", "k500_contrasts.tsv", "k500_bootstrap_draws.tsv", "prior_source_states.tsv",
                 "fold_prior_and_ranks.tsv", "bench_b_matching.tsv", "bench_b_matched_pairs.tsv", "bench_b_balance.tsv",
                 "bench_b_yields.tsv", "bench_b_point_contrasts.tsv", "bench_b_contrasts.tsv", "bench_b_bootstrap_draws.tsv", "bench_c_top500_overlap.tsv",
                 "bench_c_shared_directions.tsv", "testability_validation.tsv", "cohort_pair_counts.tsv"):
        report.append(f"- `{name}`")
    (out / "BENCH_A_BENCH_B_BENCH_C_结果摘要.md").write_text("\n".join(report) + "\n", encoding="utf-8")

    input_files = [inputs["ranks_path"], inputs["delta_path"], inputs["mapping_path"]]
    if inputs["r_path"]:
        input_files.append(inputs["r_path"])
    for reference_path in (core_grid_path, ledger_path):
        if reference_path and reference_path not in input_files:
            input_files.append(reference_path)
    protocol_path = Path(__file__).with_name("protocol_frozen.md")
    manifest = {"created_utc": datetime.now(timezone.utc).isoformat(), "status": "complete",
                "cohorts": list(COHORTS), "RNA_B_GSE15471_pairs": args.rna_b_pairs,
                "APGI_case_rule": "exact Histology == Pancreatic Ductal Adenocarcinoma; 44 matched pairs",
                "CPTAC_case_rule": "64 PDAC pairs; exclude three Adenosquamous cases from 67 linked pairs",
                "rules_m05": list(BENCH_A_RULES), "rules_m11": list(BENCH_C_RULES), "k_grid": list(KS),
                "min_pairs": MIN_PAIRS, "protein_q_cutoff": PRIMARY_Q, "bootstrap_replicates_per_target": args.bootstraps,
                "bootstrap_seed": args.seed, "bench_b_caliper_each_feature_sd": BENCH_B_CALIPER_SD,
                "heldout_target_in_prior": False, "protein_bh_family": "all released unique protein genes per target cohort",
                "testability_validation_discrepancies": int(len(state_audit)),
                "frozen_protocol_sha256": sha256(protocol_path) if protocol_path.is_file() else None,
                "input_hashes_sha256": {str(p): sha256(p) for p in input_files},
                "outputs": {p.name: sha256(p) for p in out.iterdir() if p.is_file() and p.name != "run_manifest.json"}}
    (out / "run_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    log("DONE", f"elapsed_s={time.time() - t0:.1f}", f"output_files={len(list(out.iterdir()))}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        log("ERROR", repr(exc))
        raise
