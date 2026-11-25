# Frozen methods for the benchmark analyses

This file records the fixed benchmark rules used by `benchmark_rule_contrasts.py`.
The primary RNA sets are recalculated from the paired RNA cohorts after the
documented GSE15471 quality-control pair is omitted. Historical RNA sets may be
supplied solely for a membership comparison.

The strict APGI cohort contains the 44 matched cases with exact histology
`Pancreatic Ductal Adenocarcinoma`. Strict CPTAC uses the 64 exact-label PDAC
pairs among 67 clinically linked pairs. RNA ranking B uses 35 GSE15471 pairs.

## Benchmark A: held-out protein testability prior

For each of the four target protein cohorts in turn, construct a gene-level
prior from the other three cohorts only. A source contributes 1 only when the
gene maps to a released protein row and has at least 10 finite paired
tumor-minus-normal deltas; the prior is the sum, 0–3. Record separately
whether each non-target source is absent from the released matrix, mapped but
below the pair threshold, or testable. Both absent and mapped-but-untestable
contribute zero, but remain distinct in audit tables; zero is not interpreted
as biological absence. No protein differential effect, p/q value, target-cohort
matrix, or target-case field may enter prior construction.

Rank the full shared A/B RNA universe, breaking ties by ascending gene symbol.
Compare single-A, single-B, equal-weight, and continuous RNA using their
calculated rankings and directions. The availability-only score/rank is the
prior and its direction is the equal-weight RNA direction. The combined rule
is the mean of the normalized midrank percentiles for the continuous-RNA rank
and availability-prior rank (50:50), retaining the continuous-RNA direction.
Evaluate k=100, 250, 500, 1000, and 2000.

At each fold, cohort, and k, report nomination slots, mapped genes, genes
testable at >=10 complete pairs, and q<0.05 direction-supported genes. Report
ITT support per 100 slots, mapped/testable rates per 100, and conditional
support among testable genes. Unmapped and mapped-but-untestable nominations
are unsupported in the ITT endpoint. Protein BH uses the full eligible cohort
protein universe, not only nominated genes.

For k=500 contrasts, use 1,000 bootstrap replicates of paired evaluator
patients within each target protein cohort. Rankings and the leave-one-cohort-
out prior remain fixed; each replicate recomputes protein testability and
paired tests/BH over the cohort universe. The primary estimand is the fixed-list
ITT supported-per-100 difference; also report conditional support differences
with their denominators. Compare availability-only with each RNA baseline,
RNA-plus-availability with each RNA baseline, and RNA-plus-availability with
availability-only. A common-evaluator-testable comparison is diagnostic and
post hoc; it does not replace ITT.

## Benchmark B: RNA repeatability versus single-A magnitude/precision

Use the recalculated R set and shared RNA universe. Define the candidate pool
as single-A eligible genes with finite d_A/SE_A; controls are outside R. Match
each R gene 1:1 without replacement to a control with the same sign of d_A.
Matching features are log(abs(d_A)) and log(SE_A), standardized over the
eligible shared RNA universe. Require each feature difference to be within
0.2 pooled SD and choose the smallest standardized Euclidean distance; exact
ties use ascending gene symbol. Visit R genes in ascending gene symbol. Protein
outcomes and target-cohort testability are excluded from matching.

Compare protein support for R genes in their equal-weight RNA direction and
matched controls in their single-A direction. Report slots, mapped, testable,
and q<0.05 same-direction support, with ITT and conditional yields. Bootstrap
1,000 paired target-cohort patients while keeping RNA eligibility, matching,
and directions fixed. Repeat within leave-one-cohort-out prior strata (0–3).
A comparison restricted to genes both testable in the held-out matrix is a
data-aware diagnostic only.

## Benchmark C: rank overlap and direction changes

For each pair among single-A, single-B, equal-weight, continuous,
double-significance, availability-only, and RNA-plus-availability rules,
report top-500 intersection, Jaccard, replacements, and same/opposite RNA
direction counts among shared genes. Explicitly include double-significance
versus continuous. Availability-only uses the equal-weight RNA direction;
RNA-plus-availability uses the continuous RNA direction. Protein outcomes do
not determine RNA direction.
