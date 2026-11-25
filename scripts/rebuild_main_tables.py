"""Format Tables 1–4 from a full run and optionally compare them to the submission."""

import argparse
import csv
import json
from pathlib import Path


COHORTS = {"Fudan": "Fudan", "CPTAC": "CPTAC-PDAC", "APGI": "APGI", "Ruijin": "Ruijin"}
PLATFORMS = {"Fudan": "Label-free MS", "CPTAC": "TMT log₂ MD abundance",
             "APGI": "DIA/MaxLFQ", "Ruijin": "Released proteome"}
RULES = {"single_A": "Single A", "single_B": "Single B", "equal_weight": "Equal weight",
         "continuous": "Continuous", "double_sig": "Double-significant-first"}
HELDOUT = {"continuous": "Continuous RNA", "availability_only": "Testability prior",
           "RNA_plus_availability": "RNA + prior"}


def lookup(path, fields, include=lambda row: True):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        rows = [row for row in csv.DictReader(stream, delimiter="\t") if include(row)]
    indexed = {tuple(row[field] for field in fields): row for row in rows}
    if len(rows) != len(indexed):
        raise ValueError(f"Duplicate {'/'.join(fields)} in {path}")
    return indexed


def pct(value):
    return f"{100 * float(value):.1f}%"


def signed(value):
    return f"{float(value):+.1f}".replace("-", "−")


def build(run):
    core, bench = run / "core", run / "benchmark"
    metadata = json.loads((core / "run_metadata.json").read_text(encoding="utf-8"))
    support = lookup(core / "four_cohort_support.tsv", ["cohort", "set"])
    uncertainty = lookup(core / "protein_patient_bootstrap_support.tsv", ["cohort", "set", "rate_name"])
    budget = lookup(core / "fixed_budget_grid.tsv", ["cohort", "rule"], lambda row: row["k"] == "500")
    heldout = lookup(bench / "point_yields.tsv", ["target_cohort", "rule"], lambda row: row["k"] == "500")
    contrasts = lookup(bench / "k500_contrasts.tsv", ["target_cohort", "contrast"])
    table1 = [["Cohort", "Role", "Pairs", "Quantification", "Harmonised features",
               "R mapped / testable", "R coverage"]]
    for letter, series, platform in (("A", "GSE28735", "GPL6244"), ("B", "GSE15471", "GPL570")):
        table1.append([series, f"RNA {letter}", str(metadata[f"RNA_{letter}_pairs"]),
                       f"{platform} log₂", f"{metadata[f'RNA_{letter}_genes']:,}", "—", "—"])
    for cohort, label in COHORTS.items():
        summary = metadata["protein_cohorts"][cohort]
        row = support[cohort, "R_protocol"]
        counts = f"{int(row['mapped']):,} / {int(row['testable_ge10']):,}"
        table1.append([label, "Protein", str(summary["pairs"]), PLATFORMS[cohort],
                       f"{summary['proteins']:,}", counts, pct(row["testable_rate"])])

    table2 = [["Cohort", "RNA set", "Supported / testable", "Conditional (95% CI)",
               "Unconditional (95% CI)", "Coverage (95% CI)"]]
    for name, label in (("R_protocol", "R"), ("R_strict", "R_strict")):
        for cohort, display in COHORTS.items():
            row = support[cohort, name]
            numerator, testable, total = (int(row[k]) for k in ("supported_q05", "testable_ge10", "total_frozen_genes"))
            metrics = [("supported_per_testable", numerator / testable),
                       ("supported_per_total", numerator / total),
                       ("testable_per_total", testable / total)]
            cells = []
            for metric, value in metrics:
                limits = uncertainty[cohort, name, metric]
                low, high = 100 * float(limits["lo95"]), 100 * float(limits["hi95"])
                cells.append(f"{pct(value)} ({low:.1f} to {high:.1f})")
            table2.append([display, label, f"{numerator:,} / {testable:,}", *cells])

    table3 = [["RNA ranking rule", *COHORTS.values()]]
    for rule, label in RULES.items():
        result = [budget[cohort, rule] for cohort in COHORTS]
        table3.append([label, *(f"{float(row['supported_per100']):.1f} ({row['supported']})" for row in result)])

    table4 = [["Held-out matrix", "Rule", "Supported / testable", "Conditional", "Unconditional",
               "Coverage", "ITT gain vs continuous (95% CI)"]]
    for cohort, display in COHORTS.items():
        reference = int(heldout[cohort, "continuous"]["supported_q05"])
        for rule, label in HELDOUT.items():
            row = heldout[cohort, rule]
            n, testable = int(row["supported_q05"]), int(row["testable_ge10"])
            shift = "Reference"
            if rule != "continuous":
                limits = contrasts[cohort, f"{rule}_minus_continuous"]
                low, high = signed(limits["ITT_lo95"]), signed(limits["ITT_hi95"])
                shift = f"{signed((n - reference) / 5)} ({low} to {high})"
            table4.append([display, label, f"{n} / {testable}", pct(n / testable) if testable else "—",
                           pct(n / 500), pct(testable / 500), shift])
    return [table1, table2, table3, table4]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, default=Path("build"))
    parser.add_argument("--out", type=Path, default=Path("build/main_tables"))
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    references = Path(__file__).resolve().parents[1] / "results/main_tables"
    for number, rows in enumerate(build(args.run), 1):
        name = f"Table{number}.tsv"
        with (args.out / name).open("w", encoding="utf-8", newline="") as stream:
            csv.writer(stream, delimiter="\t", lineterminator="\n").writerows(rows)
        if args.verify:
            with (references / name).open(encoding="utf-8-sig", newline="") as stream:
                expected = list(csv.reader(stream, delimiter="\t"))
            if expected != rows:
                differences = [i for i, pair in enumerate(zip(expected, rows), 1) if pair[0] != pair[1]]
                raise ValueError(f"{name}: mismatching row numbers {differences}; row counts {len(expected)}, {len(rows)}")
            print(f"{name}: all rows match")


if __name__ == "__main__":
    main()
