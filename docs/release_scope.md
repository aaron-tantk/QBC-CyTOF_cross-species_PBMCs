# v1.1.0: scope and analysis stages

This resource documents the saved analyses described in manuscript v6.
Dataset version 1.1.0 is published at
https://doi.org/10.5281/zenodo.23676187; the original dataset DOI identifies
version 1.0.0, not these additions. Existing legacy repository materials remain
in place. See [v1.1.0 notes](releases/v1.1.0.md) for the original code snapshot,
documentation-only finalization, and the owner's GitHub tag/release handoff.

## Distinguish datasets and stages

| Stage | Saved design or count described in v6 |
|---|---|
| Experimental QBC cohort | Three human and three rhesus PBMC donors, unstimulated |
| Separate stimulation reference | Published whole-blood CyTOF cohort: 46 human and 15 rhesus donors; Basal, CD40L and GM-CSF |
| Basal CyTOF display | 36,000 cells per species, 12,000 per donor; 10 human and 9 rhesus annotated populations |
| QBC display | 40,000 cells per species: 10,000 per donor from Stanford plus 10,000 Davis cells from H1 (human) or M3 (rhesus); these are not fourth donors |
| Basal QBC matching export | 30,000 cells per species; not the complete 40,000-cell display |
| Basal fused resource | 29,412 human and 28,578 rhesus computational pairs |
| Primary stimulation fusion | PBMC-restricted reference; 20,862 human and 14,414 rhesus fused rows |
| Whole-blood-reference sensitivity | Separate alternative matching design; 20,550 human and 13,501 rhesus fused rows |

Basal standalone clustering uses resolution 0.3. The recorded basal fused
resolution is 0.7, with 18 human/13 rhesus graph clusters summarized as 14/11
annotated populations. The stimulation primary resolution-0.5 clustering
has 16 human/11 rhesus clusters. These are different resources; a historical
resolution-named metadata column is not a substitute for the final annotation map.
This addition does not regenerate or certify every legacy basal clustering panel.

Each fused row has 174 assay measurements (134 QBC and 40 CyTOF), with 18 shared
anchors. These are computational pairs across assays and individuals, not 174
unique proteins measured on the same physical cell. The 59-marker QBC panel used
for fused clustering is a feature subset, not the full assay. The 140-row QBC
reagent table includes 134 antigen reagents and six isotype controls.

## Population labels and numerical summaries

- Figures 4 and 6 group on broad source-reference labels; Figure 5 uses fine
  source-reference labels. Cluster plurality is the Table S5 annotation audit,
  not a reassignment of per-cell source labels.
- Figure 4 uses pre-joint-CLR, CyTOF-side arcsinh values and paired reference-donor
  means. Eligibility is at least three fused rows per condition per donor and
  at least five paired donors. BH adjustment covers the tests actually eligible,
  not a fixed rectangular test universe. Saved permutation results are retained.
- Figures 5 and 6 use joint-CLR values. Figure 5 reuses four separate 174-feature
  graphs; its nine-marker score defines roots, not the graph feature set.
  Score-rooting consistency is not independent biological validation.
- Figure 6 summarizes 36 saved bootstrap-support tables into 12 population
  rewiring summaries. The 15 intracellular variables include total IκBα as well
  as 14 phospho epitopes. Majority-oriented edges are descriptive support
  summaries, not a guaranteed DAG or causal signaling map.
- Reused reference cells and fused rows are not independent donors. Figure 5
  first averages repeated links within reference cell; Figure 4 retains its
  saved fused-row donor aggregation. These methods are not silently harmonized.

## Source alignment and documentation finalization

The source snapshot `f0050c36f9f24ccad25b23c52b259c31346df92c`, pinned by
Zenodo, synchronized documentation with v6, including source labels,
stage-specific counts, intracellular-node terminology and primary/sensitivity
CLI selection. Its Python/R edits relative to the preceding source package
were explanatory docstrings/comments. This subsequent finalization changes
only documentation, citation metadata and related checksums; Python/R source
and compact result-table bytes remain identical to that pinned snapshot.
Current file hashes are in the repository's curated `SHA256SUMS`; historical
source hashes, `docs/upload_verification.json` and scientific validation
receipts describe earlier checks, not new scientific runs.

No matching, clustering, UMAP, trajectory fitting, bootstrap fitting or permutation
testing was rerun for this finalization. Large processed matrices, matching maps
and saved graph fits are supplied as six independently extractable derived-data
archives in published dataset v1.1.0. Extract all six into the same repository
root to reconstruct one `data/`; see the [download instructions](../README.md#download-and-verify-saved-data). The
available portable sources and deterministic summary routes improve coverage,
but do not establish complete end-to-end reproduction of every manuscript panel.
See [reproduction coverage](reproduction.md) and [environment limits](environment.md).
