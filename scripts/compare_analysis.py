"""Compare recomputed TSVs to the candidate's machine-readable reference tables."""

import argparse
import csv
import itertools
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SORT_KEYS = {
    "TableS12_core_input_manifest.tsv": ("file",),
    "TableS9_ranking_reconciliation.tsv": ("target_cohort", "k"),
}


def same_value(expected, actual):
    if expected == actual:
        return True
    if expected.lower() in ("", "nan") and actual.lower() in ("", "nan"):
        return True
    try:
        return math.isclose(float(expected), float(actual), rel_tol=5e-8, abs_tol=1e-10)
    except ValueError:
        return False


def compare(reference, generated):
    with reference.open(encoding="utf-8-sig", newline="") as left, generated.open(
        encoding="utf-8-sig", newline=""
    ) as right:
        expected = csv.DictReader(left, delimiter="\t")
        actual = csv.DictReader(right, delimiter="\t")
        if expected.fieldnames != [name for name in actual.fieldnames if name != "core_reference"]:
            raise ValueError(f"Headers differ in {reference.name}")
        keys = SORT_KEYS.get(reference.name)
        if keys:
            expected = sort_rows(expected, keys)
            actual = sort_rows(actual, keys)
        rows = 0
        for rows, (a, b) in enumerate(itertools.zip_longest(expected, actual), 1):
            if a is None or b is None:
                raise ValueError(f"Row count differs for {reference.name} at row {rows}")
            for key, value in a.items():
                if not same_value(value, b[key]):
                    raise ValueError(f"{reference.name} row {rows}, {key}: {value!r} vs {b[key]!r}")
    return rows


def sort_rows(rows, keys):
    normalized = []
    for row in rows:
        if "file" in keys:
            row["file"] = row["file"].replace("\\", "/")
        normalized.append(row)
    return sorted(normalized, key=lambda row: tuple(row[key] for key in keys))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, default=Path("build"))
    parser.add_argument("--stage", choices=["core", "benchmark", "supplement", "all"], default="all")
    args = parser.parse_args()
    failed = []
    with (ROOT / "docs/result_file_map.tsv").open(encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if args.stage != "all" and row["stage"] != args.stage:
                continue
            reference = ROOT / "results/supplementary_tables" / row["reference"]
            generated = args.run / row["generated"]
            try:
                print(f"{row['reference']}: {compare(reference, generated)} rows match")
            except (OSError, ValueError) as exc:
                failed.append(row["reference"])
                print(f"DIFFERENCE: {exc}")
    if failed:
        raise SystemExit(f"{len(failed)} comparisons did not match.")


if __name__ == "__main__":
    main()
