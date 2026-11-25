# Reproduction map

Run `analysis/run_core_benchmark.py` to produce `build/core`,
`build/benchmark` and `build/audits`. The full 1,000-draw bootstrap is the
default. The current four strict PDAC cohorts use paired cases as specified in
the publication methods, and GSE15471 contributes 35 RNA pairs after the
documented quality-control exclusion.

| Publication output | How to regenerate or inspect |
| --- | --- |
| Tables 1–4 | `python scripts/rebuild_main_tables.py --run build --verify`. The formatter uses core cohort metadata, support counts and fixed-list bootstrap, the budget grid, and benchmark point yields and paired intervals. |
| Figures 1–4 and S1–S6 | `python figures/make_figures.py --output-dir build/figures`. Each panel reads its aggregate TSV in `figures/source_data`. Reference PDF copies are in `results/reference_figures`. |
| Core and benchmark supplements | `docs/result_file_map.tsv` lists every reference TSV and the expected output path. `scripts/compare_analysis.py` compares reference and new runs. It allows numerical round-off and normalizes empty versus NaN. |
| Literature wording | `supplement/literature_coding.py` reads the stored strict labels and contextual classifications without article text. The public version of the broad coding table retains bibliographic IDs and decisions. |
| Tissue preparation, PRM and HPA | `supplement/run_protein_sensitivities.py` consumes the data-source manifest inputs and the recomputed RNA lists. The patient/sample inputs remain outside the repository. |
| GEO matching and approximate C-index | `supplement/geo_overlap.py` followed by `supplement/approximate_g3_score.py`. The reference aggregate outputs shipped with the release are under `results/audits`. |
| BioTExt case provenance | `supplement/check_pdc_case_links.py` uses live PDC/GDC and the locally generated CPTAC histology map. The included aggregate result records the analysed version. |
| Historical Figure S3 calculation | `supplement/historical_cptac_purity.py` reads the older processed proteome and phenotype files with the pinned HGNC reference and the historical RNA gene membership. |

Main-table rates distinguish conditional support (supported/testable),
unconditional support (supported/nominated), and coverage (testable/nominated).
Uncertainty and scale fields in the source tables retain their original
definitions. The supplement includes comparisons with other candidate rules and
historical assays. Their inputs and inference targets vary as described in the
source and method notes.
