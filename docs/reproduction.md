# Reproduction coverage and label provenance

This resource provides saved inputs/results and portable analysis sources. It is **not a claim of exact end-to-end reproduction of every manuscript statement**. No scientific results were recomputed with changed labels or methods. Saved data are published at [dataset DOI 10.5281/zenodo.23676187](https://doi.org/10.5281/zenodo.23676187); follow the [six-part download and extraction instructions](../README.md#download-and-verify-saved-data) before running the commands below. Validation results described here are historical saved checks, not reruns performed during GitHub documentation finalization.

## Coverage map

| Component | Available materials and status |
|---|---|
| QBC workflow/technical demonstration | Raw files remain in the original archive where provided. No reconstruction from sequencing reads or claim of reproducing every technical-demonstration panel is added here. |
| Modality-specific clustering and frequencies | Existing `Seurat-CyTOF-clustering.Rmd`, `Seurat-QBC-clustering.Rmd`, marker tables and legacy processed objects remain unchanged. These contain historical paths/settings; full agreement with final manuscript counts and normalization is not certified. Manuscript v6 distinguishes the basal CyTOF display (12,000 cells per donor; 36,000 per species), QBC display and QBC matching export; see the [stage inventory](release_scope.md). |
| Basal cross-modal fusion | Original MARIO notebook and basal fused R Markdown remain as legacy entrypoints. The new stimulation exports must not be substituted for the unstimulated Figure 3 resource. |
| Stimulation input preparation/matching | Portable recovered sources, saved matching indices/embeddings and corrected primary/sensitivity fused exports are provided. Matching was not rerun for this update. Exact refits require the original environment and random state. |
| Match-quality/comparator validation | Saved cluster, anchor-correlation and cross-method concordance tables are included. Availability of each recovered comparator source is identified in the upstream instructions; no new benchmark is claimed. |
| Figure 4 response statistics | All 375 primary and 480 sensitivity paired-test effect/p/BH values verified against saved tables. They use source-reference broad lineage, pre-joint-CLR CyTOF-side values and donor means. Saved 10,000-permutation results were retained, not rerun. No reference-cell deduplication is retroactively imposed. |
| Figure 5 | All eight numerical tables reproduced from exported primary joint-CLR assays; independent validation checked scores, source values, roots, donor eligibility/effects, rankings and all 134 marker associations per species. Saved graph fits are reused. |
| Figure 6 | All 12 population rewiring summaries reproduced exactly from 36 archived bootstrap-support tables, including majority-direction deduplication. Original network bootstraps are not rerun. |

## Source-label provenance and plurality audit

Consistent with manuscript v6, `source_reference_population` and `source_reference_broad_lineage` preserve annotations of the linked CyTOF reference cells. `audit_plurality_population` is obtained by assigning each fused cluster its most frequent source label. It is provided only for auditing.

In the primary objects, fine source and plurality labels differ in 4,283/20,862 human rows and 6,781/14,414 rhesus rows. There are eight human source populations but seven plurality populations. Changing to plurality labels would alter population memberships and may change scientific results. The source labels and saved results are retained; the audit is not a downstream label-reassignment step.

Figure 5 uses fine source labels (including `Classical_monocytes`, `Activated_dendritic_cells` and `Activated_T_cells`). Figures 4 and 6 use the source broad-lineage grouping (e.g. `Monocytes`, `DC`, `T_activated`). These label vocabularies must not be conflated.

## Assay and matching conventions

- Every released matrix starts with `cell_id`, followed by 174 measurements: 134 `-QBC` and 40 `-CyTOF`.
- `*_pre_joint_CLR.csv.gz` contains the values before the additional joint CLR step. These are **not raw counts/intensities**: QBC input is already per-donor CLR transformed, and CyTOF input is arcsinh transformed with cofactor 5.
- `*_joint_CLR.csv.gz` contains the per-row joint CLR assay used for Figure 5 and Figure 6. For a fused row x, the checked operation is `log1p(x / exp(sum(log1p(x))/174))`.
- Metadata align exactly by `cell_id`. `reference_id` identifies reused CyTOF reference cells. Donors use neutral codes with no identity lookup in these derived exports.
- `qbc_row_id`/`cytof_row_id` are zero-based row positions in the unshuffled source tables. Matching-map `row_idx_0`/`row_idx_1` refer instead to the shuffled MARIO inputs. Do not join these columns directly without reconstructing the source shuffle/filter.
- Primary and whole-blood-sensitivity data must not be pooled. The same QBC cell can appear across the two alternative designs.
- `matching/` includes the original match-index and embedding tables. Embedding_0 is QBC-side and embedding_1 CyTOF-side, already aligned in matched-pair row order; index reconstruction is specified in the assembly source.

## Figure 5

The score averages species-standardized values for nine nonconstant markers: CyTOF Ki67, HLA-DR, pS6, p4EBP1 and pERK1/2; QBC HLA-DR, CD71, CD25 and CD69. `HLA-DR-QBC` is the actual QBC feature. pSTAT5 is not a score component.

Only Basal and GM-CSF enter score scaling, fitting and summaries. Repeated links are averaged within unique reference cell before donor/population/condition averaging. A donor needs at least three unique reference cells per condition; population ranking needs at least five paired donors. Effects are descriptive.

Four independently fitted graphs use all 174 joint-CLR measurements, 50 PCs and the saved Monocle3/igraph implementation. They are **not** joint-condition nine-feature graphs. High-score top-5% roots anchor early pseudotime. Finite coverage: human Basal 5,250/7,740; human GM-CSF 5,920/6,958; rhesus Basal 4,507/4,807; rhesus GM-CSF 4,646/4,989. Independent graphs do not supply a common quantitative time axis.

Optional refit after preparation:

```bash
python code/figure5/run_graph.py --rscript YOUR_RSCRIPT --rlib YOUR_R_LIBRARY \
  --code code/figure5/fit_trajectories.R --input work/inputs \
  --out work/graphs_refit --graph human_Basal
```

Repeat for `human_GM-CSF`, `rhesus_Basal`, `rhesus_GM-CSF`. The runner requires Unix resource-limit support, a 64 MiB stack and one BLAS/OpenMP thread. The fitter checks Monocle3 1.3.1/igraph 1.6.0 and locally substitutes seed 42 for hard-coded internal seed 2016. It does not modify installed packages. Refit outputs are distinct from archived authoritative fits.

## Figure 6

The 15 intracellular measurements comprise 14 phospho epitopes and total IκBα (`IkBa.amino.terminal`). Historical variable names containing `phospho` are retained for compatibility; they do not make IκBα a phospho-specific measurement.

The archived tables contain bootstrap strength and conditional direction support, often both orientations for one unordered pair. The summary chooses the maximum-support direction per pair; ties retain the first maximum in file order. Rewiring counts unique pSTAT5-connected pairs gained, lost or reversed, with reversals counted once.

The inference source filters strength strictly above 0.85. Four archived pairs (human B GM-CSF, human CD8 CD40L, human NK GM-CSF, rhesus DC CD40L) serialize at exactly 0.85. The summaries preserve these archived rows; they do not refilter them. Featured edge totals remain human monocytes 41→42 and rhesus DC 25→27.

The original bootstrap seed/state was not recorded in the recovered inference source. Optional refitting uses `fit_networks.R PROCESSED_PRIMARY OUTPUT`; it is not an exact numerical reproduction guarantee. Majority-oriented support edges are not the same artifact as bnlearn's acyclic averaged network, and directions do not establish causal signaling.
