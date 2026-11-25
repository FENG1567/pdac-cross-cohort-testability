import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np


parser = argparse.ArgumentParser(description="Recompute the GEO specimen overlap audit.")
parser.add_argument("--inputs", required=True, type=Path)
parser.add_argument("--out", required=True, type=Path)
args = parser.parse_args()
OUT, BASE = args.out, args.inputs
OUT.mkdir(parents=True, exist_ok=True)


def load(name):
    path = BASE / f'{name}_series_matrix.txt.gz'
    with path.open('rb') as stream:
        digest = hashlib.file_digest(stream, 'sha256').hexdigest()
    meta = {}
    features = {}
    with gzip.open(path, 'rt', encoding='utf-8', newline='') as stream:
        for row in csv.reader(stream, delimiter='\t'):
            if not row:
                continue
            if row[0] == '!series_matrix_table_begin':
                break
            if row[0] in ('!Sample_title', '!Series_summary', '!Sample_characteristics_ch1'):
                meta.setdefault(row[0], []).append(row[1:])
        columns = next(csv.reader(stream, delimiter='\t'))[1:]
        for row in csv.reader(stream, delimiter='\t'):
            if not row:
                continue
            if row[0] == '!series_matrix_table_end':
                break
            if len(row) != len(columns) + 1 or row[0] in features:
                raise ValueError('Unexpected GEO probe ID or column structure')
            features[row[0]] = [float(v) if v not in ('NA', 'null', '') else np.nan for v in row[1:]]
    return columns, meta, features, digest


a, ma, da, hash_a = load('GSE28735')
b, mb, db, hash_b = load('GSE62452')
shared = sorted(da.keys() & db.keys())
X = np.asarray([da[p] for p in shared])
Y = np.asarray([db[p] for p in shared])
complete = np.isfinite(X).all(axis=1) & np.isfinite(Y).all(axis=1)


def normalized(z):
    z = z[complete]
    centered = z - z.mean(axis=0)
    return centered / np.linalg.norm(centered, axis=0)


correlation = normalized(X).T @ normalized(Y)
abest = correlation.argmax(axis=1)
bbest = correlation.argmax(axis=0)
assert len(a) == 90 and len(b) == 130
same_column = np.arange(len(a))
assert np.all(bbest[:len(a)] == same_column)
assert np.all(abest == same_column)
score = correlation[same_column, same_column]
assert any('GSE28735' in x and 'earlier' in x.lower() for row in mb['!Series_summary'] for x in row)


def tissue_labels(metadata):
    rows = metadata['!Sample_characteristics_ch1']
    tissue, = [row for row in rows if all(x.lower().startswith('tissue: ') for x in row)]
    labels = {'pancreatic tumor': 'tumor', 'adjacent pancreatic non-tumor': 'adjacent'}
    return [labels[x.lower().split(': ', 1)[1]] for x in tissue]


types = tissue_labels(mb)
assert len(types) == len(b)
a_titles, = ma['!Sample_title']
b_titles, = mb['!Sample_title']

with (OUT / 'GSE28735_GSE62452_sample_map_recomputed.tsv').open('w', encoding='utf-8', newline='') as stream:
    writer = csv.writer(stream, delimiter='\t')
    writer.writerow(['gse28735_col', 'gse28735_gsm', 'gse28735_title', 'gse62452_col', 'gse62452_gsm', 'gse62452_title', 'pearson_r'])
    for index in same_column:
        writer.writerow([index, a[index], a_titles[index], index, b[index], b_titles[index], f'{score[index]:.8f}'])

result = {
    'matching_samples': len(a),
    'zero_shared_accession_ids': not bool(set(a) & set(b)),
    'mutual_best_matching': True,
    'shared_probes': len(shared),
    'globally_complete_shared_probes': int(complete.sum()),
    'matched_r_median': float(np.median(score)),
    'matched_r_min': float(score.min()),
    'all_best_match_r_median': float(np.median(correlation.max(axis=0))),
    'all_best_match_r_min': float(correlation.max(axis=0).min()),
    'GSE62452_tissue_all': dict(Counter(types)),
    'GSE62452_tissue_overlapping': dict(Counter(types[:len(a)])),
    'GSE62452_tissue_independent': dict(Counter(types[len(a):])),
    'series_sha256': {'GSE28735': hash_a, 'GSE62452': hash_b},
}
(OUT / 'geo_match_summary.json').write_text(json.dumps(result, indent=2) + '\n', encoding='utf-8')
print(json.dumps(result, indent=2))
