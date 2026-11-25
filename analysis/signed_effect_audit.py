"""Summarise directional protein effects and candidate eligibility at the primary threshold."""

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def main(ledger_path, mapping_path, target):
    ledger = pd.read_csv(ledger_path, sep="\t")
    mapping = pd.read_csv(mapping_path, sep="\t")
    df = ledger.merge(mapping[["cohort", "gene", "reason"]], on=["cohort", "gene"], how="left", validate="one_to_one")
    rows = []
    eligibility = []
    for cohort, group in df.groupby("cohort", sort=False):
        lists = [("R", group.in_R, "frozen_R_direction"), ("R_strict", group.in_R_strict, "frozen_R_direction")]
        rules = ["single_A", "single_B", "equal_weight", "continuous", "double_sig"]
        lists.extend((rule, group[f"rank_{rule}"] <= 500, f"direction_{rule}") for rule in rules)
        for rule, selected, direction in lists:
            nominated = group[selected]
            measured = nominated[nominated.protein_testable]
            values = measured.protein_delta * measured[direction]
            assert np.isfinite(values).all(), (cohort, rule)
            quantiles = np.percentile(values, [25, 50, 75])
            individual = measured[(measured.protein_q < 0.05) & (values > 0)]
            rows.append(dict(cohort=cohort, list=rule, nominated=len(nominated), testable=len(measured),
                             supported=len(individual), conditional=len(individual)/len(measured),
                             unconditional=len(individual)/len(nominated), coverage=len(measured)/len(nominated),
                             mean_signed_protein_difference=values.mean(), q25_signed=quantiles[0],
                             median_signed=quantiles[1], q75_signed=quantiles[2], positive_signed=int((values > 0).sum()),
                             estimand="RNA-signed protein paired mean differences across testable genes, cohort input scale"))
            unmapped = ~nominated.protein_mapped
            insufficient = nominated.protein_mapped & (nominated.protein_n_pairs < 10)
            other = nominated.protein_mapped & ~insufficient & ~nominated.protein_testable
            counts = [int(unmapped.sum()), int(insufficient.sum()), int(other.sum()), len(measured)]
            assert sum(counts) == len(nominated), (cohort, rule, counts)
            eligibility.append(dict(cohort=cohort, list=rule, nominated=len(nominated),
                                    not_mapped_to_retained_protein_row=counts[0],
                                    mapped_below_minimum_complete_pairs=counts[1],
                                    mapped_other_non_estimable=counts[2], testable=counts[3],
                                    supported=len(individual)))
    target.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(target / "signed_protein_effects.tsv", sep="\t", index=False)
    pd.DataFrame(eligibility).to_csv(target / "nomination_eligibility.tsv", sep="\t", index=False)
    print(f"Verified {len(rows)} list and cohort combinations and decomposed nomination losses.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    for name in ("ledger", "mapping", "output"):
        parser.add_argument(f"--{name}", type=Path, required=True)
    args = parser.parse_args()
    main(args.ledger, args.mapping, args.output)
