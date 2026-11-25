"""Summarize the stored coding or run the lexical rule on locally obtained XML."""

import argparse
import csv
import json
import re
from collections import Counter
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TABLES = ROOT / "results/supplementary_tables"
AWARE = re.compile(r"sample overlap|overlap with GSE28735|not counted as an independent|shared samples", re.I)
STRICT = re.compile(r"independent", re.I)
GENERIC = re.compile(r"validat", re.I)


def load_rows(path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def xml_label(path):
    if not path.is_file():
        return "no_full_text"
    text = re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", path.read_text(encoding="utf-8")))
    snippets = []
    for token in ("GSE62452", "GSE28735"):
        snippets.extend(match[0].strip() for match in re.finditer(r"[^.]{0,320}" + token + r"[^.]{0,320}\.", text))
    snippets = list(dict.fromkeys(snippets))[:14]
    if not snippets:
        return "no_full_text"
    if AWARE.search(" ".join(snippets)):
        return "aware"
    both = [text for text in snippets if "GSE62452" in text and "GSE28735" in text]
    if any(STRICT.search(text) for text in both):
        return "explicit_independent"
    if any(GENERIC.search(text) for text in both):
        return "generic_validation_both"
    return "unclear"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--xml-dir", type=Path, help="Optional directory of PMC accession-named XMLs")
    args = parser.parse_args()
    stored = load_rows(TABLES / "TableS2_strict_coding.tsv")
    counts = Counter(row["label"] for row in stored)
    explicit = counts["explicit_independent"]
    report = {
        "sampled": len(stored),
        "analysable": len(stored) - counts["no_full_text"],
        "strict_hits": explicit,
        "broad_hits": explicit + counts["generic_validation_both"],
        "stored_labels": counts,
        "context_rechecks": len(load_rows(TABLES / "TableS2_contextual_rechecks.tsv")),
    }
    if args.xml_dir is not None:
        regular_accession = re.compile(r"PMC[0-9]+")
        changes = []
        for row in stored:
            pmcid = row["pmcid"]
            if not regular_accession.fullmatch(pmcid):
                raise ValueError(f"Unexpected PMC accession {pmcid!r}")
            alternative = xml_label(args.xml_dir / f"{pmcid}.xml")
            if alternative != row["label"]:
                changes.append({"pmcid": pmcid, "stored": row["label"], "local_xml": alternative})
        report["local_xml_label_differences"] = changes
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
