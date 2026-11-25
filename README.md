# Fixed-budget protein support for RNA nominations in PDAC

Reproduction materials for the **2026-11 analysis release**. The study ranks candidates using two paired RNA cohorts and evaluates protein support at fixed nomination counts in Fudan, CPTAC, APGI and Ruijin releases. The repository supplies the analysis code, published-candidate result tables, and plotting data.

| Path | Contents |
| --- | --- |
| `analysis/` | RNA input processing and strict cohort evaluation, fixed and re-ranked bootstrap, coverage decomposition, matched comparison and held-out testability benchmark |
| `supplement/` | BioTExt, ProteinTalks/HPA and Ruijin PRM sensitivity analyses, GEO overlap, approximate G3 score, PDC/GDC case-overlap check and literature coding |
| `figures/` | Plotting code and ten aggregate plotting tables |
| `results/main_tables/` | The four main tables from the analysis release |
| `results/supplementary_tables/` | Supplementary summaries and gene-level outputs, including stored coding with publication identifiers |
| `results/rna_reference/` | The two final RNA lists, directions and all ranks. These can be used as frozen supplementary analysis inputs. |
| `results/reference_figures/` | The ten figure PDFs used for the submission version |
| `results/audits/` | Additional aggregate provenance and sensitivity outputs |
| `docs/` | Input manifest, source and licensing notes, and correspondence between scripts and results |

## Replot the figures and inspect stored results

Use Python 3.12. Create a virtual environment if needed and install dependencies.

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
python figures/make_figures.py --output-dir build/figures
python supplement/literature_coding.py
```

Activate `.venv` before installing dependencies, or use its Python executable directly. On Unix-like systems use `.venv/bin/python`. On Windows use `.venv\Scripts\python.exe`. `requirements-reference.txt` records the reference Windows package versions. Differences in fonts, PDF metadata, dependency versions and platform can change the binary appearance of rendered files. The plotting TSVs and checked-in PDFs anchor the final candidate.

## Recompute the primary results

The computational starting points are the original *released* RNA and protein matrices. Run from the repository root. The input manifest gives expected SHA-256 digests, accession numbers, source URLs and local paths. Record the original study citations when reusing its data.

```bash
python scripts/inputs.py list --group core
python scripts/inputs.py fetch --group core
python scripts/inputs.py check --group core
python analysis/run_core_benchmark.py --input-root inputs --output-root build
python scripts/rebuild_main_tables.py --run build --verify
python scripts/compare_analysis.py --run build --stage core
python scripts/compare_analysis.py --run build --stage benchmark
```

If an upstream site changes a file or its URL, recover the checked release via the cited source and place it at the listed path. A digest mismatch is reported and prevents the downloaded file from being installed. The public CPTAC patient annotation and the two CCT headers are joined locally to generate the 64-case selection. The historical RNA comparator sets are optional audit inputs because the current RNA sets are recalculated from the two GEO cohorts. See [sources and coverage](docs/DATA_SOURCES.md).

The run generates subject-level intermediates under `build/` for paired testing and resampling. The repository ignores `build/`, `inputs/` and `runs/` by default. Keep those local when sharing a release. Direct outputs `build/core/` and `build/benchmark/` include the scientific source tables used to construct the main and supplementary results. See [result mappings](docs/RESULT_MAP.md).

## Run the supplemental checks

Acquire the specified separate inputs before running each analysis. For ProteinTalks the precise saved HTML has to be made available locally under its source terms. Supplemental reference results and their plots can be inspected from this repository immediately.

```bash
python scripts/inputs.py fetch --group supplement
python supplement/run_protein_sensitivities.py --assembly build/core --component biotext
python supplement/run_protein_sensitivities.py --assembly build/core --component prm
python supplement/run_protein_sensitivities.py --assembly build/core --component hpa
python supplement/check_pdc_case_links.py --case-map build/derived/cptac_case_histology.tsv
python scripts/inputs.py fetch --group geo
python supplement/geo_overlap.py --inputs inputs --out build/geo
python supplement/approximate_g3_score.py --inputs inputs --out build/geo --map build/geo/GSE28735_GSE62452_sample_map_recomputed.tsv
python scripts/inputs.py fetch --group purity
python supplement/historical_cptac_purity.py --inputs inputs --out build/historical_purity
python scripts/compare_analysis.py --run build --stage supplement
```

The PDC/GDC check queries the two live public services and reports only aggregate linkage results. GEO match keys and individual-level files produced during the audit are local working data. The 90 literature records, their lexical labels and the stored subset classifications can be summarized offline. A lexical re-analysis can optionally use locally obtained full-text XML; full-text availability can change and the stored manual classification remains the reference. For instructions and source constraints see [sources and coverage](docs/DATA_SOURCES.md).

When repeating a supplemental calculation independently, the path
`--assembly results/rna_reference` uses the same RNA sets and ranks as the
full reference analysis.

## Licensing and citation

Original code is [MIT](LICENSE). Author-owned results and figures are [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/). The frozen HGNC nomenclature reference retains its original CC0 terms. See [license scope](LICENSING.md) and [citation metadata](CITATION.cff).
