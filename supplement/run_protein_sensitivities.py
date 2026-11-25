#!/usr/bin/env python3
"""Refresh BioTExt, ProteinTalks/HPA, and Ruijin PRM supplements.

The script reads the frozen RNA signatures and rankings, but keeps
all publisher-side source inputs in their existing read-only locations. It
writes aggregate outputs only to the isolated revision output directory.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import t as student_t


ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT / "inputs"
DEFAULT_OUT = ROOT / "build/supplement"
WORKBOOK = INPUTS / "biotext/12014_2024_9450_MOESM2_ESM.xlsx"
DEFAULT_HPA = INPUTS / "hpa/hpa_pancreas.html"
DEFAULT_HGNC = INPUTS / "hgnc_complete_set.txt"
DEFAULT_CLINICAL = INPUTS / "node_ruijin"
SEED = 20261123
BOOTSTRAPS = 20000


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], columns: list[str] | None = None) -> None:
    if columns is None:
        columns = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def bh(pvalues):
    p = np.asarray(pvalues, dtype=float)
    q = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    if not ok.any():
        return q
    idx = np.flatnonzero(ok)
    order = idx[np.argsort(p[ok], kind="mergesort")]
    ranked = p[order] * len(order) / np.arange(1, len(order) + 1)
    q[order] = np.minimum.accumulate(ranked[::-1])[::-1].clip(0, 1)
    return q


def exact_signflip_p(values, null=0.0):
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)] - null
    n = len(x)
    if n < 2:
        return float("nan")
    observed = float(x.mean())
    if n <= 20:
        bits = np.arange(1 << n, dtype=np.uint32)[:, None]
        signs = (((bits >> np.arange(n, dtype=np.uint32)) & 1).astype(np.int8) * 2 - 1)
        permuted = signs @ x / n
        return float(np.mean(permuted >= observed - 1e-12))
    rng = np.random.default_rng(SEED + n)
    signs = rng.choice([-1.0, 1.0], size=(200000, n))
    permuted = signs @ x / n
    return float((1 + np.sum(permuted >= observed - 1e-12)) / (len(permuted) + 1))


def bootstrap_interval(values, rng):
    x = np.asarray(values, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 2:
        return float("nan"), float("nan")
    indices = rng.integers(0, len(x), size=(BOOTSTRAPS, len(x)))
    means = x[indices].mean(axis=1)
    return tuple(float(v) for v in np.quantile(means, [0.025, 0.975]))


def find_file(candidates: list[Path], label: str) -> Path:
    for path in candidates:
        if path.is_file():
            return path
    raise FileNotFoundError(f"Could not find {label}; checked: " + ", ".join(map(str, candidates)))


def clean_symbol(value) -> str:
    s = str(value or "").strip().upper()
    return "" if s in {"", "NA", "NAN", "NULL", "---", "-"} else s


def read_signature(path: Path) -> dict[str, int]:
    rows = read_tsv(path)
    if not rows:
        raise ValueError(f"Empty signature: {path}")
    dcol = "direction" if "direction" in rows[0] else "dir"
    result = {}
    for row in rows:
        gene = clean_symbol(row.get("gene"))
        try:
            direction = int(float(row.get(dcol, "nan")))
        except (TypeError, ValueError):
            continue
        if not gene:
            continue
        if direction not in (-1, 1):
            raise ValueError(f"Invalid direction {direction} for {gene} in {path}")
        if gene in result:
            raise ValueError(f"Duplicate symbol {gene} in {path}")
        result[gene] = direction
    if not result:
        raise ValueError(f"No valid gene/direction rows in {path}")
    return result


def read_top500(path: Path) -> dict[str, dict[str, int]]:
    rows = read_tsv(path)
    if not rows:
        raise ValueError(f"Empty rank table: {path}")
    columns = set(rows[0])
    by_rule: dict[str, list[dict]] = {}
    if {"rule", "gene", "rank", "direction"}.issubset(columns):
        for row in rows:
            try:
                rank = int(float(row["rank"]))
                direction = int(float(row["direction"]))
            except (TypeError, ValueError):
                continue
            if rank <= 500:
                if direction not in (-1, 1):
                    raise ValueError(f"Invalid direction in top 500 row: {row}")
                by_rule.setdefault(row["rule"].strip(), []).append(
                    {"gene": clean_symbol(row["gene"]), "direction": direction, "rank": rank}
                )
    elif "gene" in columns:
        rank_rules = sorted(col[len("rank_"):] for col in columns if col.startswith("rank_"))
        for rule in rank_rules:
            rank_col, direction_col = f"rank_{rule}", f"direction_{rule}"
            if direction_col not in columns:
                raise ValueError(f"Missing {direction_col} paired with {rank_col} in {path}")
            for row in rows:
                try:
                    rank = int(float(row[rank_col]))
                except (TypeError, ValueError):
                    continue
                if rank > 500:
                    continue
                try:
                    direction = int(float(row[direction_col]))
                except (TypeError, ValueError):
                    raise ValueError(f"Missing direction for top-500 {rule} row: {row['gene']}")
                if direction not in (-1, 1):
                    raise ValueError(f"Invalid direction in top 500 {rule} row: {row}")
                by_rule.setdefault(rule, []).append(
                    {"gene": clean_symbol(row["gene"]), "direction": direction, "rank": rank}
                )
    else:
        raise ValueError(f"Unsupported rank table schema in {path}: {list(rows[0])}")
    out = {}
    for rule, selected in sorted(by_rule.items()):
        selected.sort(key=lambda row: row["rank"])
        genes = {row["gene"]: row["direction"] for row in selected if row["gene"]}
        if len(selected) != 500 or len(genes) != 500:
            raise ValueError(f"Rule {rule!r}: expected 500 ranked genes, got {len(selected)} rows / {len(genes)} unique")
        out[rule] = genes
    if not out:
        raise ValueError(f"No rule has 500 ranks in {path}")
    return out


def read_biotext_matrix(workbook: Path):
    xls = pd.ExcelFile(workbook, engine="openpyxl")
    try:
        proteins = pd.read_excel(xls, sheet_name="Proteomics_data", engine="openpyxl")
        summary = pd.read_excel(xls, sheet_name="Summary", header=None, engine="openpyxl")
    finally:
        xls.close()
    source_desc = " ".join(str(x) for x in summary.to_numpy().ravel() if pd.notna(x))
    columns = list(proteins.columns)
    if "GeneSymbols" not in columns or len(columns) <= 4:
        raise ValueError("BioTExt publisher workbook schema does not match the expected Proteomics_data sheet")
    sample_cols = columns[4:]
    symbols = proteins["GeneSymbols"].fillna("").astype(str).str.strip()
    ambiguous = symbols.str.contains(r";|,|\\", regex=True)
    valid = symbols.ne("") & ~ambiguous
    data = proteins.loc[valid, ["GeneSymbols"] + sample_cols].copy()
    data["GeneSymbols"] = data["GeneSymbols"].map(clean_symbol)
    for col in sample_cols:
        data[col] = pd.to_numeric(data[col], errors="coerce")
    gene_area = data.groupby("GeneSymbols", sort=True)[sample_cols].sum(min_count=1)
    logmat = np.log2(gene_area.where(gene_area > 0))
    return logmat, sample_cols, source_desc, len(proteins), len(columns), int(valid.sum())


def paired_deltas(logmat: pd.DataFrame, tissue: str):
    sample_cols = list(logmat.columns)
    cases = sorted({str(col).split("_")[0] for col in sample_cols if "_" in str(col)})
    pairs = []
    for case in cases:
        tumor, nat = f"{case}_Tumor_{tissue}", f"{case}_NAT_{tissue}"
        if tumor in logmat.columns and nat in logmat.columns:
            pairs.append((case, tumor, nat))
    genes = logmat.index
    diffs = pd.DataFrame(index=genes, columns=[p[0] for p in pairs], dtype=float)
    for case, tumor, nat in pairs:
        diffs[case] = logmat[tumor] - logmat[nat]
    return pairs, diffs


def run_biotext(out: Path, workbook: Path, signature_paths: dict[str, Path], ranks_path: Path):
    logmat, _, source_desc, protein_rows, protein_cols, valid_rows = read_biotext_matrix(workbook)
    sets = {name: read_signature(path) for name, path in signature_paths.items()}
    ranked_sets = read_top500(ranks_path)
    for rule, genes in ranked_sets.items():
        sets[f"heldout_{rule}"] = genes
    tissues = ["Epithelial", "Bulk", "Stroma"]
    summary_rows = []
    per_group = {}
    rng = np.random.default_rng(SEED)
    for list_name, direction_map in sets.items():
        genes = list(direction_map)
        directions = pd.Series(direction_map, dtype=float)
        for tissue in tissues:
            pairs, tissue_diffs = paired_deltas(logmat, tissue)
            available = [gene for gene in genes if gene in tissue_diffs.index]
            if available:
                deltas = tissue_diffs.loc[available].reindex(genes)
            else:
                deltas = pd.DataFrame(index=genes, columns=tissue_diffs.columns, dtype=float)
            tested_mask = deltas.notna().sum(axis=1) >= min(3, len(pairs)) if pairs else pd.Series(False, index=genes)
            tested_genes = list(tested_mask[tested_mask].index)
            signed = deltas.loc[tested_genes].mul(directions.loc[tested_genes], axis=0) if tested_genes else pd.DataFrame()
            pscore = signed.mean(axis=0, skipna=True) if len(signed) else pd.Series(dtype=float)
            evaluable = signed.notna().sum(axis=0) if len(signed) else pd.Series(dtype=float)
            support = signed.gt(0).sum(axis=0) / evaluable.replace(0, np.nan) if len(signed) else pd.Series(dtype=float)
            scores = pscore.dropna().to_numpy(dtype=float)
            support_values = support.dropna().to_numpy(dtype=float)
            lo, hi = bootstrap_interval(scores, rng)
            formal = list_name in {"R_protocol", "R_strict"}
            p_mean = exact_signflip_p(scores) if formal else float("nan")
            p_support = exact_signflip_p(support_values, 0.5) if formal else float("nan")
            if tested_genes:
                medians = deltas.loc[tested_genes].median(axis=1, skipna=True)
                signed_medians = medians * directions.loc[tested_genes]
                concordant = int((signed_medians > 0).sum())
                discordant = int((signed_medians < 0).sum())
                tied = int((signed_medians == 0).sum())
                gene_fraction = concordant / (concordant + discordant) if concordant + discordant else float("nan")
            else:
                concordant = discordant = tied = 0
                gene_fraction = float("nan")
            summary_rows.append({
                "list": list_name,
                "tissue": tissue,
                "source_list_size": len(genes),
                "matched_unambiguous_gene_symbols": len(available),
                "tested_genes_min3_pairs": len(tested_genes),
                "gene_coverage_pct": 100 * len(tested_genes) / len(genes) if genes else float("nan"),
                "matched_case_pairs": len(pairs),
                "median_evaluable_genes_per_pair": float(np.median(evaluable.to_numpy())) if len(evaluable) else float("nan"),
                "mean_patient_signed_log2FC": float(np.mean(scores)) if len(scores) else float("nan"),
                "bootstrap95_low": lo,
                "bootstrap95_high": hi,
                "mean_patient_directional_support_fraction": float(np.mean(support_values)) if len(support_values) else float("nan"),
                "gene_median_direction_concordance_fraction": gene_fraction,
                "gene_concordant": concordant,
                "gene_discordant": discordant,
                "gene_tied": tied,
                "p_exact_signflip_signed_mean_gt0": p_mean,
                "p_exact_signflip_support_gt0_5": p_support,
            })
            per_group[(list_name, tissue)] = {"deltas": deltas, "tested_genes": tested_genes, "pairs": pairs, "directions": directions}
    primary = [row for row in summary_rows if row["list"] in {"R_protocol", "R_strict"}]
    qvalues = bh([row["p_exact_signflip_signed_mean_gt0"] for row in primary])
    qmap = {(row["list"], row["tissue"]): float(q) for row, q in zip(primary, qvalues)}
    for row in summary_rows:
        row["q_BH_primary6"] = qmap.get((row["list"], row["tissue"]), float("nan"))

    comparison_rows = []
    for list_name in sets:
        ep = per_group[(list_name, "Epithelial")]
        bulk = per_group[(list_name, "Bulk")]
        shared = sorted(set(ep["tested_genes"]) & set(bulk["tested_genes"]))
        epi_scores, bulk_scores, epi_support, bulk_support = [], [], [], []
        shared_cases = sorted(set(p[0] for p in ep["pairs"]) & set(p[0] for p in bulk["pairs"]))
        direction_map = sets[list_name]
        for case in shared_cases:
            if not shared:
                continue
            e = ep["deltas"].loc[shared, case] * pd.Series(direction_map).loc[shared]
            b = bulk["deltas"].loc[shared, case] * pd.Series(direction_map).loc[shared]
            ok = e.notna() & b.notna()
            if ok.any():
                epi_scores.append(float(e[ok].mean()))
                bulk_scores.append(float(b[ok].mean()))
                epi_support.append(float((e[ok] > 0).mean()))
                bulk_support.append(float((b[ok] > 0).mean()))
        delta = np.asarray(epi_scores) - np.asarray(bulk_scores)
        delta_support = np.asarray(epi_support) - np.asarray(bulk_support)
        lo, hi = bootstrap_interval(delta, rng)
        support_lo, support_hi = bootstrap_interval(delta_support, rng)
        comparison_rows.append({
            "list": list_name,
            "contrast": "Tumor_Epithelial_vs_Tumor_Bulk_on_shared_genes",
            "shared_tested_genes": len(shared),
            "paired_cases": len(delta),
            "mean_patient_score_difference_Epithelial_minus_Bulk": float(delta.mean()) if len(delta) else float("nan"),
            "bootstrap95_low": lo,
            "bootstrap95_high": hi,
            "p_exact_signflip_difference_gt0": exact_signflip_p(delta),
            "mean_patient_support_fraction_difference": float(delta_support.mean()) if len(delta_support) else float("nan"),
            "support_fraction_bootstrap95_low": support_lo,
            "support_fraction_bootstrap95_high": support_hi,
        })

    # Gene-wise paired t tests use all testable protein genes in each tissue for BH.
    gene_test_rows = []
    for tissue in tissues:
        pairs, tissue_diffs = paired_deltas(logmat, tissue)
        records = []
        for gene in tissue_diffs.index:
            values = tissue_diffs.loc[gene].dropna().to_numpy(dtype=float)
            if len(values) < 10:
                continue
            mean = float(values.mean())
            centered = values - mean
            ss = float(np.dot(centered, centered))
            if ss == 0:
                pvalue = 1.0 if mean == 0 else 0.0
            else:
                statistic = abs(mean) / math.sqrt(ss / (len(values) * (len(values) - 1)))
                pvalue = float(2 * student_t.sf(statistic, len(values) - 1))
            records.append({"gene": gene, "n_pairs": len(values), "mean_log2FC_tumor_minus_NAT": mean, "p_value": pvalue})
        qvalues = bh([row["p_value"] for row in records])
        for row, q in zip(records, qvalues):
            row["q_BH"] = float(q)
        test_lookup = {row["gene"]: row for row in records}
        for set_name in ("R_protocol", "R_strict"):
            direction_map = sets[set_name]
            subset = [test_lookup[g] for g in direction_map if g in test_lookup]
            concordant = sum(row["q_BH"] < 0.05 and row["mean_log2FC_tumor_minus_NAT"] * direction_map[row["gene"]] > 0 for row in subset)
            genes_n = len(direction_map)
            evaluable_n = len(subset)
            gene_test_rows.append({
                "set": set_name,
                "compartment": tissue,
                "matched_pairs": len(pairs),
                "BH_family_size": len(records),
                "genes": genes_n,
                "matrix_evaluable": evaluable_n,
                "concordant_BH05": int(concordant),
                "conditional_pct": 100 * concordant / evaluable_n if evaluable_n else float("nan"),
                "unconditional_pct": 100 * concordant / genes_n if genes_n else float("nan"),
                "imputed_matrix_evaluable_pct": 100 * evaluable_n / genes_n if genes_n else float("nan"),
            })

    write_tsv(out / "TableS7_biotext_directional_summary.tsv", summary_rows)
    write_tsv(out / "TableS7_biotext_compartment_comparison.tsv", comparison_rows)
    write_tsv(out / "TableS7_biotext_gene_tests.tsv", gene_test_rows)
    return {
        "publisher_workbook": str(workbook),
        "source_description": source_desc,
        "source_declares_imputation": "imput" in source_desc.lower(),
        "proteomics_rows": protein_rows,
        "proteomics_columns": protein_cols,
        "unambiguous_gene_symbol_rows": valid_rows,
        "unique_unambiguous_symbols": len(logmat.index),
        "compartments": {t: {"matched_pairs": len(paired_deltas(logmat, t)[0])} for t in tissues},
        "ranked_rule_counts": {rule: len(genes) for rule, genes in ranked_sets.items()},
        "summary_rows": summary_rows,
        "gene_test_rows": gene_test_rows,
    }


def hpa_support_limits(point, direction):
    cutoff = -math.log10(0.05)
    display_half_unit = 0.0005
    x = point["x"] * direction
    y = point["y"]
    definitely = x - display_half_unit >= 0.5 and y - display_half_unit > cutoff
    possibly = x + display_half_unit >= 0.5 and y + display_half_unit > cutoff
    return definitely, possibly


def run_hpa(out: Path, hpa_page: Path, hgnc_path: Path, signature_paths: dict[str, Path], ranks_path: Path):
    page = hpa_page.read_text(encoding="utf-8")
    plot_key = "$('#cancer_guo_volcano_Pancreatic_carcinoma').volcanoPlot("
    if page.count(plot_key) != 1:
        raise ValueError("Expected one ProteinTalks pancreatic carcinoma volcano plot")
    points, _ = json.JSONDecoder().raw_decode(page.split(plot_key, 1)[1])
    by_id, symbol_rows = {}, {}
    for point in points:
        match = re.fullmatch(r"/(ENSG\d+)-([^/]+)/cancer", point["u"])
        if not match or match[2] != point["g"]:
            raise ValueError(f"Malformed ProteinTalks point: {point.get('u')}")
        if match[1] in by_id:
            raise ValueError(f"Duplicate HPA Ensembl identifier: {match[1]}")
        if not math.isfinite(point["x"]) or not math.isfinite(point["y"]):
            raise ValueError("Non-finite HPA point")
        if point["ml"] != 0:
            raise ValueError("Unexpected HPA volcano label setting")
        by_id[match[1]] = point
        symbol_rows.setdefault(point["g"].upper(), []).append(point)
    for symbol, rows in symbol_rows.items():
        if len({(row["x"], row["y"], row["p"]) for row in rows}) != 1:
            raise ValueError(f"Repeated HPA symbol has non-identical endpoints: {symbol}")
    by_symbol = {symbol: rows[0] for symbol, rows in symbol_rows.items()}

    approved = {}
    for row in read_tsv(hgnc_path):
        ensembl = row.get("ensembl_gene_id", "")
        if row.get("status") == "Approved" and re.fullmatch(r"ENSG\d+", ensembl):
            symbol = clean_symbol(row.get("symbol"))
            if symbol in approved:
                raise ValueError(f"Duplicate approved HGNC symbol: {symbol}")
            approved[symbol] = ensembl
    if len(approved) < 18000:
        raise ValueError(f"HGNC crosswalk unexpectedly small: {len(approved)}")

    def tabulate(label, nominated):
        eligible = [gene for gene in nominated if approved.get(gene) in by_id]
        limits = [hpa_support_limits(by_id[approved[gene]], nominated[gene]) for gene in eligible]
        low, high = sum(a for a, _ in limits), sum(b for _, b in limits)
        denominator, evaluable = len(nominated), len(eligible)
        return {
            "set_or_rule": label,
            "candidate_count": denominator,
            "HPA_tested_displayed_count": evaluable,
            "candidates_with_approved_ensembl_id": sum(gene in approved for gene in nominated),
            "symbol_linked_displayed_count_for_comparison": len(set(nominated) & set(by_symbol)),
            "supported_count_rounding_low": low,
            "supported_count_rounding_high": high,
            "conditional_pct_low": 100 * low / evaluable if evaluable else float("nan"),
            "conditional_pct_high": 100 * high / evaluable if evaluable else float("nan"),
            "unconditional_pct_low": 100 * low / denominator if denominator else float("nan"),
            "unconditional_pct_high": 100 * high / denominator if denominator else float("nan"),
            "HPA_tested_displayed_coverage_pct": 100 * evaluable / denominator if denominator else float("nan"),
        }

    signature_rows = [tabulate(name, read_signature(path)) for name, path in signature_paths.items()]
    ranked = read_top500(ranks_path)
    budget_rows = [tabulate(rule, genes) for rule, genes in ranked.items()]
    write_tsv(out / "TableS8_hpa_support.tsv", signature_rows)
    write_tsv(out / "TableS8_hpa_fixed_budget_k500.tsv", budget_rows)
    metadata = {
        "dataset": "ProteinTalks-DB pancreatic carcinoma, paired proDA",
        "source": "https://www.proteinatlas.org/humanproteome/cancer/pancreatic+cancer",
        "source_sha256": sha256(hpa_page),
        "hgnc_crosswalk_sha256": sha256(hgnc_path),
        "rankings_sha256": sha256(ranks_path),
        "frozen_set_sha256": {name: sha256(path) for name, path in signature_paths.items()},
        "hgnc_approved_symbols_with_ensembl_id": len(approved),
        "pairs_reported_by_HPA": 37,
        "volcano_source_rows": len(points),
        "distinct_ensembl_ids": len(by_id),
        "distinct_symbols": len(by_symbol),
        "duplicates_identical_on_displayed_effect_q": len(points) - len(by_symbol),
        "bounds": "Display precision only; not confidence intervals.",
        "limitations": [
            "Published paired proDA summary only; patient-level observations were unavailable.",
            "The published minimum detection is 10% per group; detection coverage cannot be independently recomputed.",
            "Case-level PDAC histology and cohort overlap cannot be independently adjudicated from this plot.",
            "Candidate symbols were linked via approved HGNC IDs to plotted Ensembl identifiers.",
        ],
    }
    (out / "hpa_public_input_qc.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    return {"metadata": metadata, "support": signature_rows, "fixed_budget": budget_rows}


def sample_key(value):
    label = str(value or "").strip().upper().replace("_", "-")
    hit = re.search(r"RJ-(\d+)-(\d+)(?=$|[^0-9])", label)
    return f"RJ-{int(hit[1])}-{int(hit[2])}" if hit else label


def read_xlsx_rows(path: Path):
    book = pd.ExcelFile(path, engine="openpyxl")
    try:
        frame = pd.read_excel(book, sheet_name=book.sheet_names[0], header=None, engine="openpyxl")
    finally:
        book.close()
    return frame.where(pd.notna(frame), None).values.tolist()


def run_prm(out: Path, clinical_dir: Path, signature_paths: dict[str, Path], ranks_path: Path):
    prm_file = clinical_dir / "20230403_PDAC_PRM_exp.xlsx"
    meta_file = clinical_dir / "20230403_PDAC_PRM_meta.xlsx"
    dia_file = clinical_dir / "20221227_PDAC_PRO_exp.xlsx"
    prm, meta = read_xlsx_rows(prm_file), read_xlsx_rows(meta_file)
    dia_columns = pd.read_excel(dia_file, engine="openpyxl", header=None, nrows=1).iloc[0].tolist()
    if not prm or not meta:
        raise ValueError("Ruijin PRM or metadata workbook is empty")
    prm_keys = [sample_key(value) for value in prm[0][1:] if value]
    dia_keys = {sample_key(value) for value in dia_columns[1:] if value}
    meta_head = [str(value or "").strip().lower() for value in meta[0]]
    idx_id, idx_status, idx_time = meta_head.index("id"), meta_head.index("os_status"), meta_head.index("os_month")
    meta_rows = [(sample_key(row[idx_id]), row[idx_status], row[idx_time]) for row in meta[1:] if row[idx_id]]
    meta_keys = [row[0] for row in meta_rows]
    matched = [row for row in meta_rows if row[0] in set(prm_keys)]
    usable = [row for row in matched if isinstance(row[2], (float, int)) and math.isfinite(row[2]) and row[2] > 0 and row[1] is not None]
    targets = {clean_symbol(row[0]) for row in prm[1:] if row and row[0]}
    targets.discard("")
    signatures = {name: set(read_signature(path)) for name, path in signature_paths.items()}
    ranked = read_top500(ranks_path)
    sets = {**signatures, **ranked}
    rows = []
    for name, genes in sets.items():
        overlap = sorted(targets & set(genes))
        rows.append({
            "set_or_rule": name,
            "candidate_count": len(genes),
            "PRM_panel_target_count": len(targets),
            "overlap_count": len(overlap),
            "overlap_symbols": ";".join(overlap),
            "pct_PRM_panel_covered": 100 * len(overlap) / len(targets) if targets else float("nan"),
            "passes_three_target_gate": len(overlap) >= 3,
        })
    write_tsv(out / "TableS10_ruijin_prm_overlap.tsv", rows)
    report = {
        "sample_label_family": {
            name: sum(key.startswith("RJ-") for key in values)
            for name, values in (("PRM", prm_keys), ("DIA", dia_keys), ("OS", meta_keys))
        },
        "PRM_sample_columns": len(prm_keys),
        "PRM_unique_keys": len(set(prm_keys)),
        "DIA_unique_patient_keys": len(dia_keys),
        "OS_records": len(meta_rows),
        "OS_unique_keys": len(set(meta_keys)),
        "PRM_OS_direct_matches": len(set(prm_keys) & set(meta_keys)),
        "PRM_DIA_direct_matches": len(set(prm_keys) & dia_keys),
        "OS_usable_matched_records": len(usable),
        "OS_status_codes_in_usable_records": dict(sorted(Counter(str(row[1]) for row in usable).items())),
        "PRM_gene_count": len(targets),
        "PRM_gene_labels": sorted(targets),
        "overlap_sensitivity": rows,
        "refreshed_R_intersection": sorted(targets & signatures["R_protocol"]),
        "refreshed_R_strict_intersection": sorted(targets & signatures["R_strict"]),
        "rank_rule_top500_intersections": {name: sorted(targets & set(genes)) for name, genes in ranked.items()},
        "GSE15471_QC_note": "RNA signatures/ranks supplied by the refreshed assembly after the RNA-B GSE15471 QC pair exclusion.",
        "hashes": {
            prm_file.name: sha256(prm_file),
            meta_file.name: sha256(meta_file),
            dia_file.name: sha256(dia_file),
            ranks_path.name: sha256(ranks_path),
            **{path.name: sha256(path) for path in signature_paths.values()},
        },
        "interpretation": "Aggregate eligibility/overlap only. The PRM data have no tumor-adjacent paired specimens; no survival association was evaluated.",
    }
    (out / "prm_eligibility.json").write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--assembly", type=Path, default=ROOT / "build/core")
    parser.add_argument("--workbook", type=Path, default=WORKBOOK)
    parser.add_argument("--hpa-page", type=Path, default=DEFAULT_HPA)
    parser.add_argument("--hgnc", type=Path, default=DEFAULT_HGNC)
    parser.add_argument("--clinical-dir", type=Path, default=DEFAULT_CLINICAL)
    parser.add_argument("--component", choices=["all", "biotext", "hpa", "prm"], default="all")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    signature_paths = {
        name: find_file([
            args.assembly / "rna_qc_rankings" / f"{name}.tsv",
            args.assembly / f"{name}.tsv",
        ], name)
        for name in ("R_protocol", "R_strict")
    }
    ranks_path = find_file([
        args.assembly / "rna_qc_rankings" / "rna_scores_all_ranks.tsv",
        args.assembly / "rna_scores_all_ranks.tsv",
    ], "rna_scores_all_ranks.tsv")
    if args.component == "biotext":
        run_biotext(args.out, args.workbook, signature_paths, ranks_path)
        return
    if args.component == "hpa":
        run_hpa(args.out, args.hpa_page, args.hgnc, signature_paths, ranks_path)
        return
    if args.component == "prm":
        run_prm(args.out, args.clinical_dir, signature_paths, ranks_path)
        return
    source_inputs = [args.workbook, args.hpa_page, args.hgnc,
                     args.clinical_dir / "20230403_PDAC_PRM_exp.xlsx",
                     args.clinical_dir / "20230403_PDAC_PRM_meta.xlsx",
                     args.clinical_dir / "20221227_PDAC_PRO_exp.xlsx", *signature_paths.values(), ranks_path]
    missing = [str(path) for path in source_inputs if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing required supplemental input(s): " + ", ".join(missing))

    biotext = run_biotext(args.out, args.workbook, signature_paths, ranks_path)
    hpa = run_hpa(args.out, args.hpa_page, args.hgnc, signature_paths, ranks_path)
    prm = run_prm(args.out, args.clinical_dir, signature_paths, ranks_path)
    outputs = [
        "TableS7_biotext_directional_summary.tsv",
        "TableS7_biotext_gene_tests.tsv",
        "TableS7_biotext_compartment_comparison.tsv",
        "TableS8_hpa_support.tsv",
        "TableS8_hpa_fixed_budget_k500.tsv",
        "TableS10_ruijin_prm_overlap.tsv",
        "hpa_public_input_qc.json",
        "prm_eligibility.json",
    ]
    manifest = {
        "run": "Refreshed BioTExt directional/BH gene tests, HPA ProteinTalks summary, and Ruijin PRM overlap",
        "run_date": "2026-11-23",
        "seed": SEED,
        "bootstraps": BOOTSTRAPS,
        "RNA_B_GSE15471_QC_pair": "excluded upstream in the refreshed assembly; signatures/ranks are consumed as frozen inputs",
        "method": {
            "BioTExt_directional": "RNA-direction-weighted paired log2 tumor-minus-NAT effects; at least 3 pairs per gene; one-sided exact patient-pair sign-flip tests for R_protocol/R_strict; BH across the six formal list-by-compartment tests; top-500 rule lists are descriptive.",
            "BioTExt_gene_tests": "Two-sided paired one-sample t tests on gene-summed log2 intensity differences; at least 10 complete pairs; BH over all testable protein genes within each compartment; concordance counted after BH.",
            "HPA": "Published paired proDA volcano summary with rounding-aware display bounds; approved HGNC-to-Ensembl linkage.",
            "Ruijin_PRM": "Aggregate fixed-panel target overlap and sample-linkage QC only; no prognostic association.",
        },
        "inputs_sha256": {path.name: sha256(path) for path in source_inputs},
        "outputs_sha256": {name: sha256(args.out / name) for name in outputs},
        "biotext": {k: v for k, v in biotext.items() if k not in {"summary_rows", "gene_test_rows"}},
        "hpa": hpa,
        "prm_summary": {k: v for k, v in prm.items() if k not in {"PRM_gene_labels", "overlap_sensitivity", "refreshed_R_intersection", "refreshed_R_strict_intersection", "rank_rule_top500_intersections"}},
    }
    (args.out / "supplemental_refresh_manifest.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "out": str(args.out),
        "tables": outputs[:6],
        "biotext_directional_rows": len(biotext["summary_rows"]),
        "biotext_gene_test_rows": len(biotext["gene_test_rows"]),
        "hpa_support": hpa["support"],
        "hpa_fixed_budget": hpa["fixed_budget"],
        "prm_overlap": prm["overlap_sensitivity"],
        "manifest": str(args.out / "supplemental_refresh_manifest.json"),
    }, indent=2, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
