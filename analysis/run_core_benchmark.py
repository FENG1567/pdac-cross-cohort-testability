#!/usr/bin/env python3
"""Run the frozen core analysis, audits, and benchmark module.

Expected input-root contents (historical RNA sets are optional):
  GSE28735_series_matrix.txt.gz, GSE15471_series_matrix.txt.gz,
  GPL6244.annot.gz, GPL570.annot.gz,
  13045_2022_1384_MOESM25_ESM.xlsx,
  cptac_orig/proteomics_gene_level_MD_abundance_{tumor,normal}.cct,
  apgi/protein_matrix.csv, hgnc_complete_set.txt,
  node_ruijin/20221227_PDAC_PRO_exp.xlsx.

The CPTAC public clinical file is supplied with --cptac-clinical (or found
under one of the documented candidate names in input-root). The pipeline
derives its pair-level histology mask at run time; no case rows are bundled.
Requires numpy, pandas, scipy, and openpyxl.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from pathlib import Path


HERE = Path(__file__).resolve().parent
DEFAULT_REPO_ROOT = HERE.parent


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--input-root", required=True, type=Path,
                        help="local directory containing the public/staged inputs listed above")
    parser.add_argument("--output-root", required=True, type=Path,
                        help="new or dedicated output directory; core, rankings, audits, and benchmark are created below it")
    parser.add_argument("--repo-root", type=Path, default=DEFAULT_REPO_ROOT,
                        help=f"package root used as the subprocess working directory (default: {DEFAULT_REPO_ROOT})")
    parser.add_argument("--cptac-clinical", type=Path,
                        help="original public CPTAC clinical JSON/TSV; defaults to input-root/cptac_orig/data_clinical_patient.txt when present")
    parser.add_argument("--cptac-mapping-output", type=Path,
                        help="derived case mapping path (default: output-root/derived/cptac_case_histology.tsv)")
    parser.add_argument("--historical-sets", type=Path,
                        help="optional historical comparison directory containing R_protocol.tsv and R_strict.tsv; primary RNA sets are recalculated either way")
    parser.add_argument("--bootstraps", type=int, default=1000,
                        help="benchmark bootstrap replicates per target cohort (default: 1000)")
    parser.add_argument("--seed", type=int, default=20261123,
                        help="benchmark random seed (default: 20261123)")
    parser.add_argument("--skip-audits", action="store_true",
                        help="skip the standalone matrix-scale and signed-effect aggregate audits")
    return parser.parse_args()


def resolve(path: Path) -> Path:
    return path.expanduser().resolve()


def inside(child: Path, parent: Path) -> bool:
    try:
        child.relative_to(parent)
        return True
    except ValueError:
        return False


def discover_clinical(input_root: Path) -> Path:
    candidates = (
        "cptac_orig/data_clinical_patient.txt",
        "cptac_clinical.json", "cptac_clinical.tsv", "clinical_attributes.json",
        "clinical_attributes.tsv", "clinical.json", "clinical.tsv",
        "cptac_orig/clinical.json", "cptac_orig/clinical.tsv",
        "cptac_orig/clinical_attributes.json", "cptac_orig/clinical_attributes.tsv",
    )
    for relative in candidates:
        path = input_root / relative
        if path.is_file():
            return path.resolve()
    raise FileNotFoundError("could not auto-detect public CPTAC clinical JSON/TSV; pass --cptac-clinical")


def require_files(paths: list[Path]) -> None:
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError("required input file(s) missing:\n  " + "\n  ".join(missing))


def run_step(label: str, command: list[str], cwd: Path) -> None:
    print(f"[{label}] {json.dumps(command, ensure_ascii=False)}", flush=True)
    subprocess.run(command, cwd=cwd, check=True)


def main() -> None:
    args = parse_args()
    repo_root = resolve(args.repo_root)
    input_root = resolve(args.input_root)
    output_root = resolve(args.output_root)
    if not input_root.is_dir():
        raise FileNotFoundError(f"input-root is not a directory: {input_root}")
    if not repo_root.is_dir():
        raise FileNotFoundError(f"repo-root/package-root is not a directory: {repo_root}")
    if input_root == output_root or inside(input_root, output_root) or inside(output_root, input_root):
        raise ValueError("input-root and output-root must be separate, non-overlapping directories")
    if args.bootstraps < 1:
        raise ValueError("--bootstraps must be at least 1")
    if not (HERE / "core_final.py").is_file() or not (HERE / "benchmark_rule_contrasts.py").is_file():
        raise FileNotFoundError("portable core/benchmark scripts must remain beside this entry point")

    clinical = resolve(args.cptac_clinical) if args.cptac_clinical else discover_clinical(input_root)
    tumor = input_root / "cptac_orig/proteomics_gene_level_MD_abundance_tumor.cct"
    normal = input_root / "cptac_orig/proteomics_gene_level_MD_abundance_normal.cct"
    mapping = (resolve(args.cptac_mapping_output) if args.cptac_mapping_output
               else output_root / "derived" / "cptac_case_histology.tsv")
    if inside(mapping, input_root):
        raise ValueError("derived CPTAC mapping must be written outside input-root")
    require_files([clinical, tumor, normal])

    required_input = [
        input_root / "GSE28735_series_matrix.txt.gz",
        input_root / "GSE15471_series_matrix.txt.gz",
        input_root / "GPL6244.annot.gz",
        input_root / "GPL570.annot.gz",
        input_root / "13045_2022_1384_MOESM25_ESM.xlsx",
        input_root / "apgi" / "protein_matrix.csv",
        input_root / "hgnc_complete_set.txt",
        input_root / "node_ruijin" / "20221227_PDAC_PRO_exp.xlsx",
    ]
    require_files(required_input)
    # Historical RNA sets are an optional comparator and are never inferred from
    # a local step2 directory; the primary sets are recalculated from RNA data.
    historical_sets = resolve(args.historical_sets) if args.historical_sets else None
    if historical_sets is not None:
        require_files([historical_sets / "R_protocol.tsv", historical_sets / "R_strict.tsv"])

    core_out = output_root / "core"
    audit_out = output_root / "audits"
    benchmark_out = output_root / "benchmark"
    output_root.mkdir(parents=True, exist_ok=True)
    audit_out.mkdir(parents=True, exist_ok=True)

    run_step("derive-cptac-case-map", [sys.executable, str(HERE / "derive_cptac_mapping.py"),
        "--clinical", str(clinical), "--tumor", str(tumor), "--normal", str(normal),
        "--output", str(mapping)], repo_root)

    core_command = [sys.executable, str(HERE / "core_final.py"),
        "--inputs", str(input_root), "--out", str(core_out), "--cptac-histology", str(mapping)]
    if historical_sets is not None:
        core_command.extend(["--historical-sets", str(historical_sets)])
    run_step("core", core_command, repo_root)

    if not args.skip_audits:
        run_step("matrix-scale-audit", [sys.executable, str(HERE / "matrix_scale_audit.py"),
            "--inputs", str(input_root), "--output", str(audit_out / "matrix_scale_audit.tsv")], repo_root)
        run_step("signed-effect-audit", [sys.executable, str(HERE / "signed_effect_audit.py"),
            "--ledger", str(core_out / "analysis_ledger.tsv"),
            "--mapping", str(core_out / "protein_mapping_testability.tsv"),
            "--output", str(audit_out)], repo_root)

    run_step("benchmark", [sys.executable, str(HERE / "benchmark_rule_contrasts.py"),
        "--core-dir", str(core_out), "--out-dir", str(benchmark_out),
        "--rna-b-pairs", "35", "--bootstraps", str(args.bootstraps), "--seed", str(args.seed)], repo_root)
    print(f"Pipeline complete: core={core_out}; audits={audit_out}; benchmark={benchmark_out}; "
          f"derived_cptac_mapping={mapping}", flush=True)


if __name__ == "__main__":
    main()
