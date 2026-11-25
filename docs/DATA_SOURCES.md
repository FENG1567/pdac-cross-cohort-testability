# Data sources and reuse

The analysis starts with processed public releases. File names, repository
accessions, URLs, local destinations and reference digests are in
[`input_sources.tsv`](input_sources.tsv). Source data remain under their
depositors' terms and should be cited when reused.

| Input | Source and version notes |
| --- | --- |
| RNA nominations | NCBI GEO [GSE28735](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE28735) and [GSE15471](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE15471), with GPL6244 and GPL570. The GSE15471 exclusion is reconciled with [KEGGdzPathwaysGEO](https://bioconductor.org/packages/KEGGdzPathwaysGEO). The public GEO array identifiers used for that QC rule are part of the implementation; the repo does not bundle an expression matrix. |
| Fudan | Tong et al., [doi:10.1186/s13045-022-01384-3](https://doi.org/10.1186/s13045-022-01384-3), Supplementary Table S3A, workbook `13045_2022_1384_MOESM25_ESM.xlsx`. Also indexed at [Figshare](https://doi.org/10.6084/m9.figshare.21628888.v1). |
| CPTAC | Cao et al., [doi:10.1016/j.cell.2021.08.023](https://doi.org/10.1016/j.cell.2021.08.023) and [PDC000270](https://pdc.cancer.gov/pdc/study/PDC000270). Use the median-normalized tumor and normal `.cct` files in the [study-specific LinkedOmics release](https://www.linkedomics.org/data_download/CPTAC-PDAC/). Use the [cBioPortal DataHub patient annotation](https://github.com/cBioPortal/datahub/tree/master/public/paad_cptac_2021) for pathology. The derived case file stays local. |
| APGI | Aref et al., [doi:10.1158/2767-9764.CRC-25-0229](https://doi.org/10.1158/2767-9764.CRC-25-0229), [PRIDE PXD059074](https://www.ebi.ac.uk/pride/archive/projects/PXD059074). Use the archive file `protein_matrix.csv` with the pinned digest. The transformation and precise measurement scale upstream of this release remain incompletely documented. |
| Ruijin | [doi:10.1038/s41591-023-02790-x](https://doi.org/10.1038/s41591-023-02790-x), NODE [OEP004100](https://www.biosino.org/node/project/detail/OEP004100). The exact three public OED records are in the manifest. The main analysis reads the processed abundance release; the PRM sensitivity uses the original panel and phenotype workbooks locally. The DIA column list follows the exact workbook header. |
| HGNC | The [HGNC complete set](https://hgnc.genenames.org/download/) (CC0) frozen at the recorded digest is bundled separately under `third_party/hgnc/`. This reference contains nomenclature and gene identifiers. |
| BioTExt | Savage et al., [doi:10.1186/s12014-024-09450-3](https://doi.org/10.1186/s12014-024-09450-3), [PDC000504](https://pdc.cancer.gov/pdc/study/PDC000504). The journal workbook is used for the imputed proteome sensitivity. PDC and GDC public case metadata support the 15-versus-64 aggregate overlap check; the script does not treat specimen- or aliquot-level independence as established. |
| ProteinTalks/HPA | [Paired-proDA pancreatic cancer display](https://www.proteinatlas.org/humanproteome/cancer/pancreatic+cancer), captured 2026-11-22, SHA-256 in the manifest. The [HPA license page](https://www.proteinatlas.org/about/licence) lists separate use terms for ProteinTalks-DB data. The reference HTML is held outside the release. For exact replay place an authorized, hash-matched copy in `inputs/hpa/hpa_pancreas.html`. The public repository provides the paper's aggregated counts and plotting data for inspection. A newly fetched HTML page may differ from the historical input. |
| GEO overlap and approximate G3 score | [GSE62452](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62452) and GSE28735. The submitters describe the shared samples. The approximate published score uses 15 distinct coefficient entries from Chen et al., [doi:10.1093/g3journal/jkab296](https://doi.org/10.1093/g3journal/jkab296), evaluated using Harrell C. Individual scores, survival values and the specimen map stay local. |
| Historical CPTAC purity | [LinkedOmicsKB Release 1, PDAC](https://kb.linkedomics.org/download), two files listed under the `purity` group. This is the historical within-tumor result shown in Figure S3. The historical gene membership table contains the RNA rule membership used by that archive for this calculation. The current fixed-budget results use the recomputed RNA sets. |
| Literature audit | [Europe PMC](https://europepmc.org/) search `GSE28735 AND GSE62452`. The stored frame contained 321 records, with 314 PMC IDs. The public TSV retains the sampled 90 record IDs and coding. If needed, acquire available PMC XMLs by ID from Europe PMC or NCBI in accordance with source terms. `supplement/literature_coding.py --xml-dir <folder>` expects filenames such as `PMC10054004.xml` and reports differences from the stored lexical labels. The contextual review was performed manually. |

The raw spectra, original matrix and clinical workbooks, third-party article
excerpts, saved ProteinTalks HTML and external case-level files remain with
their original sources or in local project storage. The repository uses
author-prepared source summaries, cohort- and gene-level results, and figures.

`results/supplementary_tables/TableS12_core_input_manifest.tsv` is the historical
ten-file measurement and annotation list. The expanded `input_sources.tsv`
additionally pins the clinical file used to derive the CPTAC histology mask
and the external inputs for the supplemental analyses.
