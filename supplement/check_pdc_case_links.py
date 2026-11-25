"""Compare open PDC and GDC case identifiers, reporting only aggregate counts."""

import argparse
import csv
import json
from pathlib import Path
from uuid import UUID
from urllib.error import HTTPError
from urllib.request import Request, urlopen


PRIMARY = "PDC000270"
BIOTEXT = "PDC000504"



def api(url, payload):
    request = Request(url, json.dumps(payload).encode(),
                      {"Content-Type": "application/json", "User-Agent": "Academic provenance check"})
    try:
        with urlopen(request, timeout=20) as response:
            result = json.load(response)
    except HTTPError as error:
        try:
            details = json.load(error)
        except ValueError:
            raise RuntimeError(f"PDC API HTTP {error.code}") from None
        messages = [entry.get("message", "") for entry in details.get("errors", [])]
        raise RuntimeError(f"API HTTP {error.code}: {messages}") from None
    if result.get("errors"):
        raise ValueError(result["errors"][0]["message"])
    return result["data"]


def query(text):
    return api("https://pdc.cancer.gov/graphql", {"query": text})


def cohort(study):
    cases = {}
    total = None
    for offset in range(0, 300, 100):
        selection = f'paginatedCasesSamplesAliquots(pdc_study_id:"{study}", offset:{offset}, limit:100)'
        data = query(f'{{ {selection} {{ total casesSamplesAliquots {{ case_id case_submitter_id samples {{ gdc_sample_id }} }} }} }}')
        rows = data["paginatedCasesSamplesAliquots"]
        total = int(rows["total"])
        cases.update({x["case_submitter_id"]: x for x in rows["casesSamplesAliquots"]})
        if len(cases) == total:
            break
    assert len(cases) == total and all(row["case_id"] for row in cases.values())
    print(study, "complete PDC case count", len(cases), flush=True)
    return cases


def sample_uuid(value):
    try:
        return str(UUID(str(value)))
    except ValueError:
        return None


def gdc_links(cases):
    samples = {key: {sample_uuid(s["gdc_sample_id"]) for s in row["samples"] or []} - {None}
               for key, row in cases.items()}
    identifiers = sorted(set().union(*samples.values())) if samples else []
    print("PDC cases with a GDC sample UUID", sum(map(bool, samples.values())), "of", len(cases), flush=True)
    if not identifiers:
        data = api("https://api.gdc.cancer.gov/cases", {"filters": {"op": "in", "content": {
            "field": "submitter_id", "value": sorted(cases)}},
            "fields": "case_id,submitter_id,project.project_id", "size": 250, "format": "JSON"})
        assert len(data["hits"]) == data["pagination"]["total"]
        lookup = {key: {row["case_id"] for row in data["hits"] if row["submitter_id"] == key} for key in cases}
        print("GDC submitter matches", sum(len(ids) == 1 for ids in lookup.values()), "of", len(cases),
              "ambiguous", sum(len(ids) > 1 for ids in lookup.values()), flush=True)
        return {key: next(iter(ids)) for key, ids in lookup.items() if len(ids) == 1}
    data = api("https://api.gdc.cancer.gov/cases", {"filters": {"op": "in", "content": {
        "field": "samples.sample_id", "value": identifiers}},
        "fields": "case_id,samples.sample_id", "size": 1000, "format": "JSON"})
    assert len(data["hits"]) == data["pagination"]["total"]
    samples_to_cases = {sample["sample_id"]: case["case_id"] for case in data["hits"]
                        for sample in case.get("samples", [])}
    links = {key: {samples_to_cases[s] for s in keys if s in samples_to_cases} for key, keys in samples.items()}
    assert all(len(ids) <= 1 for ids in links.values())
    return {key: next(iter(ids)) for key, ids in links.items() if ids}


def pdc_ids(cases):
    return {row["case_id"] for row in cases.values()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case-map", required=True, type=Path)
    args = parser.parse_args()
    primary = cohort(PRIMARY)
    bio = cohort(BIOTEXT)
    with args.case_map.open(encoding="utf-8", newline="") as stream:
        exact64 = {r["case_key"] for r in csv.DictReader(stream, delimiter="\t") if r["keep_pdac"] == "1"}
    assert len(exact64) == 64
    target = {key: primary[key] for key in sorted(exact64 & primary.keys())}
    print(json.dumps({"matched_primary_analysis_cases": len(target), "primary_analysis_cases": len(exact64),
                      "shared_pdc_case_uuid_with_main": len(pdc_ids(target) & pdc_ids(bio)),
                      "shared_pdc_case_uuid_with_release": len(pdc_ids(primary) & pdc_ids(bio))}), flush=True)
    bt_gdc = gdc_links(bio)
    main_gdc = gdc_links(target)
    print(json.dumps({"biotext_cases_with_gdc_id": len(bt_gdc), "main_cases_with_gdc_id": len(main_gdc),
                      "shared_gdc_case_uuid_with_main": len(set(bt_gdc.values()) & set(main_gdc.values()))}), flush=True)


if __name__ == "__main__":
    main()
