"""List, validate, and optionally retrieve inputs from their original sources."""

import argparse
import csv
import gzip
import hashlib
import shutil
import sys
from pathlib import Path
from urllib.error import URLError
from urllib.request import Request, urlopen


ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def matches(path, sha256):
    return path.is_file() and (not sha256 or digest(path) == sha256)


def install(source_opener, target, sha256):
    target.parent.mkdir(parents=True, exist_ok=True)
    temp = target.with_name(target.name + ".partial")
    try:
        with source_opener() as source, temp.open("wb") as output:
            shutil.copyfileobj(source, output)
        if not matches(temp, sha256):
            raise ValueError(f"SHA-256 {digest(temp)} does not match {sha256}")
        temp.replace(target)
    finally:
        temp.unlink(missing_ok=True)


def fetch(row, target):
    if row["path"] == "hgnc_complete_set.txt":
        archive = ROOT / "third_party/hgnc/hgnc_complete_set.txt.gz"
        install(lambda: gzip.open(archive, "rb"), target, row["sha256"])
    else:
        request = Request(row["direct_url"], headers={"User-Agent": "PDAC-reproduction/2026.11.24"})
        install(lambda: urlopen(request, timeout=90), target, row["sha256"])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=["list", "check", "fetch"])
    parser.add_argument("--group", choices=["core", "supplement", "geo", "purity", "all"], default="core")
    parser.add_argument("--inputs", type=Path, default=ROOT / "inputs")
    args = parser.parse_args()
    root = args.inputs.resolve()
    errors = 0
    with (ROOT / "docs/input_sources.tsv").open(encoding="utf-8", newline="") as stream:
        rows = [row for row in csv.DictReader(stream, delimiter="\t")
                if args.group == "all" or args.group in row["group"].split(",")]
    for row in rows:
        target = (root / row["path"]).resolve()
        if not target.is_relative_to(root):
            raise ValueError(f"Input destination is outside the input root: {row['path']}")
        sha256 = row["sha256"]
        if args.action == "list":
            print(row["path"], row["accession"], row["source"], sep="\t")
        elif matches(target, sha256):
            print(f"{'VERIFIED' if sha256 else 'PRESENT, unpinned'} {row['path']}")
        elif target.exists():
            errors += 1
            print(f"VERSION MISMATCH {row['path']}: {digest(target)}; expected {sha256}", file=sys.stderr)
        elif args.action == "fetch" and (row["direct_url"] or row["path"] == "hgnc_complete_set.txt"):
            try:
                fetch(row, target)
                print(f"{'VERIFIED' if sha256 else 'DOWNLOADED, unpinned'} {row['path']}")
            except (OSError, URLError, ValueError) as exc:
                errors += 1
                print(f"FETCH FAILED {row['path']}: {exc}; see {row['source']}", file=sys.stderr)
        else:
            errors += 1
            print(f"MISSING {row['path']}; see {row['source']}", file=sys.stderr)
    return int(errors > 0)


if __name__ == "__main__":
    sys.exit(main())
