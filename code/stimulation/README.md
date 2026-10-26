# Portable QBC upstream analysis source addition

This directory adds source code, **not a new scientific analysis or regenerated
release**. No remote changes, matching, clustering, UMAP, integration, model
training, or Figure 4 permutation runs were performed to prepare it. Manuscripts,
artwork, editorial reports, donor lookups, and large data are intentionally absent.

## Scientific conventions that must not be conflated

* `CyTOF_annotation` is the **saved source reference-cell fine label**.
  `broad_lineage` is its **saved source broad label**, named
  `CyTOF_broad_lineage` in fused tables. Assembly transfers both by the original
  balanced-reference row index, not by a new QBC cluster annotation.
* Figure 4 groups on those saved **broad per-cell source labels**. These labels
  are retained even where a fused cluster contains multiple source populations.
* Table S5 is a **plurality/composition audit** of the saved fine labels within a
  fused QBC cluster. Its historical column `transferred_CyTOF_label` means the
  most frequent partner label **in that audit only**. It is never written back
  to the per-cell fine or broad label. Manuscript v6 describes the source-label
  rule used by the saved downstream analyses; plurality remains an audit only.
* `annotate_cytof.R` intentionally replaces the original script's fixed
  cluster-ID-to-label reassignment with import of **saved source labels**, with
  row/cell-ID checks. New stochastic cluster IDs cannot reconstruct those labels.
  The saved annotations are required inputs, not inferred by this addition.
* Existing explicit preprocessing seeds (42) and the Figure 4 permutation seed
  (42) are preserved. **The original complete RNG state for MARIO, Seurat,
  Harmony, and CellFuse is unknown.** No new global seed was added. Library
  defaults and unseeded operations can make reruns differ from saved artifacts.
* Package versions in `provenance/dependencies_observed.json` are **observed
  installed versions during this port**, not proven original analysis versions
  and not a tested environment lock. The original complete environment is
  unavailable. Even an observed dependency commit does not establish the commit
  used for every original result.

## Layout and input contract

`scripts/` contains portable Python/R entry points and shared helpers;
`markers/` contains the unchanged, byte-for-byte Table S4 CSV (59 unique marker
names). It is a marker/clone/species-reactivity panel, not a patient table.
The original Table S4 name mapping remains:
`FcεRIα -> FcRI`, `HLA-A,B,C -> HLA-ABC`, `Integrin β7 -> integrinB7`.
Clustering retains the source intersection behavior rather than silently
substituting or renaming other markers.

Choose a project root and provide:

```
project/
  raw/cytof/Combined_Data_part1of6.csv ... Combined_Data_part6of6.csv
  raw/qbc/human_01.csv ... rhesus_03.csv
  config/qbc_inputs.csv
  data/qbc_{human,rhesus}_30k.csv
  data/cytof_{human,rhesus}_balanced.csv
  data/cytof_{human,rhesus}_balanced_annotations.csv
  mario_out/{human,rhesus}_matches.csv
  mario_out/pbmc_{human,rhesus}_matches.csv       # primary PBMC-reference analysis
```

You may start with saved prepared inputs and saved matches instead of rebuilding
raw inputs. Preserve **all row orders**: match indices are zero-based indices
after the exact source shuffle/truncation, while fused `cytof_row_id` refers to
the original balanced table and its row-aligned annotation table.

The annotation CSV must contain `donor`, `stimulation`, `CyTOF_annotation`,
`broad_lineage`; the optional `cell` column is checked when attaching to an RDS.
Full balanced expression and annotation tables must have identical row-wise
donor/stimulation metadata. Match indices must be finite, integral and in range.
Do not filter/reorder prepared inputs independently of their saved matches.
Smoke matching results are not inputs to the production assembler.

### Neutral donor labels and filenames

The QBC manifest example uses **only neutral indexed labels and paths**. Paths
are relative to the manifest unless absolute. Arrange rows in the original
published notebook's input-file order to preserve donor-ordered QBC cell rows.
Stage public raw archive filenames under indexed local names, or set the local
manifest paths yourself. No identity mapping needs to be published. This package
does not link public filename strings to participant identities.

QBC manifest labels must be `human_qbc_01`, etc. For existing original tables,
helpers privately replace sorted unique donor strings with
`{species}_{modality}_01`, etc.; a lookup is never logged or written. CyTOF
relabeling is **after** the original donor/condition grouping, sampling, and row
shuffle, so it cannot change the sampled cells. Use full donor universes when
converting linked tables; do not independently re-index filtered donor subsets.
Donor labels are metadata; expression, labels, test grouping membership and row
indices are not changed. Neutral indices are not a demographic crosswalk.

## Commands

Run from this addition's directory. `python` and `Rscript` must refer to suitable
installed environments. Set a real project directory, for example:

```bash
export QBC_PROJECT_ROOT="$PWD/project"
```

Every entry point supports `--help`. `--input-dir`/`--output-dir` override
project-root defaults. No machine-specific or home-directory input defaults are
embedded. Keep this directory's `scripts/` and `markers/` layout together.

### Saved-artifact route (no matching or clustering)

```bash
python scripts/prepare_pbmc.py
python scripts/assemble_fused.py --tag pbmc
# Or omit --tag for the full-reference analysis:
python scripts/assemble_fused.py
```

The primary PBMC-reference filter excludes exactly `Granulocyte/Myeloid`, `Platelets`,
`Erythrocytes`, `Unresolved` on the source broad labels. It preserves remaining
row order. Assembly preserves the source shuffle seed, 30,000 CyTOF truncation
for the full reference, shared-anchor within-pair correlations, and matching QC.
Output expression tables and QC are written under `project/fused/`.

The legacy no-tag default is the **whole-blood-reference sensitivity** design,
not the primary design. Use `--tag pbmc` explicitly for primary assembly,
clustering and Figure 4 statistics. Do not combine outputs from the two designs.

Optional saved-label attachment to an existing, row-aligned reference object:

```bash
Rscript scripts/annotate_cytof.R --species human \
  --input project/data/cytof_human_balanced_clustered.rds \
  --annotations project/data/cytof_human_balanced_annotations.csv \
  --output project/data/cytof_human_balanced_annotated.rds
```

### Optional future computational reruns (not executed for this addition)

These commands can be expensive and/or stochastic. Use **new output directories**
if saved release artifacts must remain unchanged.

```bash
python scripts/build_cytof_balanced.py
Rscript scripts/build_qbc_clr.R --input project/config/qbc_inputs.csv
python scripts/prepare_pbmc.py
python scripts/run_mario.py --pbmc-only
python scripts/assemble_fused.py --tag pbmc
Rscript scripts/cluster_fused.R --species human --tag pbmc
Rscript scripts/cluster_fused.R --species rhesus --tag pbmc
python scripts/fig4_stats.py --tag pbmc
```

Raw QBC preparation retains per-file Seurat CLR (`margin=2`), 134 input protein
columns and the original index-column heuristic. It does not resample raw files:
the published 10k-per-donor raw files are expected. CyTOF preparation retains
isotope mass mapping, arcsinh(x/5), donor-condition caps 260 human/800 rhesus,
source group order, and both source seeded selection/shuffle operations.

`run_mario.py` keeps every source parameter: normalization, five batches,
one match/cell, sparsities 400/400, components 10/12, five canonical correlations,
four weights, ten clusters/filter components, bad proportion 0.2,
20 filter iterations, no kNN, embedding dimension 10 and max-iteration argument
500 (50 for source `--smoke`). The installed MARIO implementation may itself
interpret embedding iteration arguments differently; no library internals were
changed. Matching has no new seed. Both modalities retain all original features,
with only the 18 named anchors overlapping.

`cluster_cytof.R` is available as an optional reference-reclustering method;
do **not** use its new cluster IDs to replace saved source labels.
`cluster_fused.R` retains CLR normalization, the Table S4 QBC-feature intersection,
20 PCs, neighbors on PCs 1:20, resolutions 0.5 and 1.2, and UMAP.
It writes RDS and Table S5 audits to `project/results/` by default (the source
put the RDS beside the fused CSV). `--rds-output-dir` can override that location.

### Figure 4 implementation, exactly as code rather than paraphrased prose

The 15 intracellular nodes are 14 phospho epitopes plus total IκBα
(`IkBa.amino.terminal`); historical `phospho` variable/output names are retained.
For each species/broad population/15-node/two-stimulus combination, donor means
are calculated on the fused CyTOF arcsinh-scale values. A donor must have at
least three cells in Basal and in that stimulus; at least five eligible donors
are needed. The paired t-test compares stimulated and basal donor means. BH
adjustment pools **all eligible tests actually run across both species**, not an
assumed complete rectangular test universe (375 eligible primary tests and
480 whole-blood-reference sensitivity tests in the saved analyses).

Cross-species correlations compare species means of eligible donor deltas:
globally across shared lineage/node pairs (at least five pairs) and separately
within a shared lineage across nodes (at least eight nodes).
Global permutations shuffle **rhesus node deltas within each lineage**; the
within-lineage tests shuffle rhesus deltas across nodes. They do not swap species
labels. The source uses 10,000 permutations and the +1 p-value correction, with
one `default_rng(42)` stream in the original iteration order. No alternative
statistics, effect definitions, labels, p-value thresholds, or seeds were added.

### Optional validation comparators

Dependency provenance and inspected APIs are recorded, but stochastic execution
was **not** verified by rerunning models. Install the indicated upstream packages
in caller-managed environments; this addition does not vendor third-party code.

Harmony (full reference only; source workflow uses the 18 shared anchors,
modality-wise scaling, 15-PC exact PCA, Harmony, CyTOF 1-NN and ARI/NMI):

```bash
# Deterministic export only, with source labels, not a new MARIO run:
python scripts/run_mario.py --prepare-only
Rscript scripts/concordance_harmony.R --species human
Rscript scripts/concordance_harmony.R --species rhesus
```

CellFuse (uses the inspected upstream `python_funcs.py`; its Python environment
is separate from the core environment):

```bash
python scripts/build_cellfuse_data.py
python scripts/build_cellfuse_pbmc.py
# Set CELLFUSE_PYTHON_DIR to an installed CellFuse checkout's python/ directory.
python scripts/run_cellfuse.py --species human --tag pbmc \
  --cellfuse-python "$CELLFUSE_PYTHON_DIR"
```

The reference train/validation split remains 70/30, stratified, random_state 42;
query anchors and source broad labels are unchanged. All source training and
prediction hyperparameters are retained, including their differences. Training
remains unseeded and date-based model naming is retained. A later model run can
produce different predictions. Comparator summary/plotting scripts without
available standalone source were not invented.

## Verification and provenance

`verification/receipt.json` records syntax, dependency-free help/missing-path
errors, Table S4 marker count, source-vs-port deterministic matching-input order,
four real 24-pair assembly comparisons (two species x full/PBMC), two narrow
Figure 4 donor-test comparisons to saved results, and real-subset CyTOF/QBC input
preparation comparisons. These checks do not regenerate scientific results.
Numerical comparisons use the explicit tolerances recorded in the receipt.
The historical source-comparison harness is not a portable release entrypoint.
Use the repository's `code/validate_release.py` to validate the packaged inputs,
matching alignment and saved Figure 5/6 numerical summaries.

`provenance/source_sha256.tsv` identifies the original source bytes and original
marker panel; the repository checksum manifest identifies delivered bytes. Hashes
prove file identity, **not** scientific reproducibility or original environment
identity. Required saved labels and matches are supplied in the candidate data
archive and are not regenerated by the source port.
