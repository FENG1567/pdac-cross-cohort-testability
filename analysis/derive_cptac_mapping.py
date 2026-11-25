#!/usr/bin/env python3
"""Derive the CPTAC paired-case histology mask from public clinical inputs.

The output is a run-time intermediate. No case-level mapping is bundled with
the analysis package.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


EXPECTED_PAIRED_CASES = 67
EXPECTED_PDAC_CASES = 64
PDAC_LABEL = "Pancreatic Ductal Adenocarcinoma"


def _records(payload):
    """Yield records from the public CPTAC/PDC JSON response shapes."""
    if isinstance(payload, list):
        yield from (row for row in payload if isinstance(row, dict))
    elif isinstance(payload, dict):
        if "clinicalAttributeId" in payload:
            yield payload
        else:
            for child in payload.values():
                if isinstance(child, (list, dict)):
                    yield from _records(child)


def _field(columns, names):
    by_normalized = {"".join(c.lower().split()): c for c in columns}
    for name in names:
        found = by_normalized.get("".join(name.lower().split()))
        if found is not None:
            return found
    return None


def read_histology(path: Path) -> dict[str, str]:
    """Read the public clinical JSON or TSV and reject duplicate case rows."""
    lookup: dict[str, str] = {}
    if path.suffix.lower() == ".json":
        payload = json.loads(path.read_text(encoding="utf-8-sig"))
        rows = [row for row in _records(payload)
                if str(row.get("clinicalAttributeId", "")).strip() == "HISTOLOGY_DIAGNOSIS"]
        for row in rows:
            case = str(row.get("patientId", "") or "").strip()
            histology = str(row.get("value", "") or "").strip()
            if not case:
                raise ValueError("public clinical JSON has a histology row without patientId")
            if case in lookup:
                raise ValueError(f"duplicate HISTOLOGY_DIAGNOSIS record for case {case}")
            lookup[case] = histology
    else:
        with path.open(encoding="utf-8-sig", newline="") as stream:
            def clinical_lines():
                for line in stream:
                    if line.startswith("#"):
                        candidate = line.lstrip("#")
                        columns = candidate.rstrip("\r\n").split("\t")
                        has_case = _field(columns, ("PATIENT_ID", "patientId", "CASE_ID", "caseId"))
                        has_histology = _field(columns, ("HISTOLOGY_DIAGNOSIS", "histologyDiagnosis"))
                        has_long_form = (_field(columns, ("clinicalAttributeId",)) and
                                         _field(columns, ("value",)))
                        if has_case and (has_histology or has_long_form):
                            yield candidate
                        continue
                    yield line

            reader = csv.DictReader(clinical_lines(), delimiter="\t")
            if not reader.fieldnames:
                raise ValueError(f"clinical TSV has no header: {path}")
            case_col = _field(reader.fieldnames, ("PATIENT_ID", "patientId", "CASE_ID", "caseId"))
            histology_col = _field(reader.fieldnames, ("HISTOLOGY_DIAGNOSIS", "histologyDiagnosis"))
            attribute_col = _field(reader.fieldnames, ("clinicalAttributeId",))
            value_col = _field(reader.fieldnames, ("value",))
            if case_col is None:
                raise ValueError("clinical TSV needs PATIENT_ID/patientId or CASE_ID/caseId")
            if histology_col is None and not (attribute_col and value_col):
                raise ValueError("clinical TSV needs HISTOLOGY_DIAGNOSIS or clinicalAttributeId/value columns")
            for row in reader:
                if attribute_col and value_col:
                    if str(row.get(attribute_col, "")).strip() != "HISTOLOGY_DIAGNOSIS":
                        continue
                    value = row.get(value_col, "")
                else:
                    value = row.get(histology_col, "")
                case = str(row.get(case_col, "") or "").strip()
                if not case:
                    raise ValueError("clinical TSV has a histology row without a case identifier")
                if case in lookup:
                    raise ValueError(f"duplicate HISTOLOGY_DIAGNOSIS record for case {case}")
                lookup[case] = str(value or "").strip()
    if not lookup:
        raise ValueError(f"no HISTOLOGY_DIAGNOSIS records found in {path}")
    return lookup


def read_sample_ids(path: Path) -> set[str]:
    """Read sample identifiers from a wide CPTAC CCT/TSV matrix header."""
    with path.open(encoding="utf-8-sig", newline="") as stream:
        for line in stream:
            fields = next(csv.reader([line.rstrip("\r\n")], delimiter="\t"), [])
            sample_ids = [value.strip() for value in fields[1:]]
            if len(sample_ids) > 20:
                if any(not value for value in sample_ids):
                    raise ValueError(f"blank sample identifier in matrix header: {path}")
                if len(set(sample_ids)) != len(sample_ids):
                    raise ValueError(f"duplicate sample identifiers in matrix header: {path}")
                return set(sample_ids)
    raise ValueError(f"no wide sample header (>20 samples) found in {path}")


def derive_mapping(clinical: Path, tumor: Path, normal: Path, output: Path) -> Counter:
    histology_by_case = read_histology(clinical)
    paired_cases = read_sample_ids(tumor) & read_sample_ids(normal)
    if len(paired_cases) != EXPECTED_PAIRED_CASES:
        raise ValueError(f"expected {EXPECTED_PAIRED_CASES} paired CPTAC cases; found {len(paired_cases)}")
    missing = paired_cases - set(histology_by_case)
    if missing:
        raise ValueError(f"public clinical join is incomplete for {len(missing)} paired cases")
    counts = Counter(histology_by_case[case] for case in paired_cases)
    pdac_count = counts.get(PDAC_LABEL, 0)
    if pdac_count != EXPECTED_PDAC_CASES:
        raise ValueError(f"expected {EXPECTED_PDAC_CASES} exact-label PDAC pairs; found {pdac_count}")

    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t", lineterminator="\n")
        writer.writerow(("case_key", "histology", "keep_pdac"))
        for case in sorted(paired_cases):
            histology = histology_by_case[case]
            writer.writerow((case, histology, int(histology == PDAC_LABEL)))
    print(f"Derived CPTAC map: paired={len(paired_cases)}; exact-label PDAC={pdac_count}; "
          f"other={len(paired_cases) - pdac_count}; output={output}", flush=True)
    return counts


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--clinical", required=True, type=Path,
                        help="original public CPTAC clinical JSON or TSV")
    parser.add_argument("--tumor", required=True, type=Path,
                        help="CPTAC tumor abundance CCT/TSV matrix")
    parser.add_argument("--normal", required=True, type=Path,
                        help="CPTAC normal abundance CCT/TSV matrix")
    parser.add_argument("--output", required=True, type=Path,
                        help="run-time derived case mapping output TSV")
    args = parser.parse_args()
    derive_mapping(args.clinical, args.tumor, args.normal, args.output)


if __name__ == "__main__":
    main()
