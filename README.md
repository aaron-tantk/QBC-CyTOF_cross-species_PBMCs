# QBC–CyTOF cross-species PBMC resource

Analysis code for human and rhesus macaque QBC/CyTOF profiling and fusion with a CyTOF cohort reference.

**Code version 1.1.0; manuscript-v6-aligned resource.** Dataset version 1.1.0 is published at [DOI 10.5281/zenodo.23676187](https://doi.org/10.5281/zenodo.23676187). This is a **dataset DOI**, not a software DOI. Cite the exact data version and code commit actually used. See [v1.1.0 notes](docs/releases/v1.1.0.md) for code provenance and the GitHub tag/release handoff.

## Study design

- QBC: PBMCs from three human and three rhesus macaque donors.
- CyTOF eference: a separate published CyTOF cohort of 46 human and 15 rhesus donors under Basal, CD40L and GM-CSF conditions.
- MARIO links cells across assays and different individuals; fused measurements are not measurements of the same physical cell. Stimulation-associated signaling transferred from the reference is not direct evidence that QBC specimens were stimulated.
- Primary stimulation fusion uses a PBMC-restricted reference: 20,862 human and 14,414 rhesus fused rows. Whole-blood-reference sensitivity data contain 20,550 and 13,501 rows.

## Label provenance and reproduction scope

Consistent with manuscript v6, the saved analyses use each matched **CyTOF source-reference cell's population label** (fine labels for Figure 5; broad labels for Figures 4 and 6). They do **not** use fused-cluster plurality labels for those population summaries. Source and plurality labels disagree for 4,283/20,862 human and 6,781/14,414 rhesus primary rows. Human plurality transfer produces seven populations, whereas the saved fine source labels produce eight.

This resource preserves the saved results rather than silently relabeling cells or recomputing statistics. Table S5 plurality labels remain available as a **cluster annotation audit**, not as the labels used for downstream population summaries. These source additions and documentation align with v6; legacy notebooks are retained, not certified as complete end-to-end reproductions. See [reproduction coverage](docs/reproduction.md) and [release scope and analysis stages](docs/release_scope.md).

## Contents

| Location | Purpose |
|---|---|
| `code/figure5/` | Two-condition score preparation, four separate 174-feature trajectories, descriptive donor summaries and marker associations |
| `code/figure6/` | Reproduce summaries of 36 archived support tables; optional inference source, not refitted for this update |
| `code/data/` | Export existing fused assay objects with explicit label provenance |
| `results/figure5/` | Compact effect and coverage summaries; all eight verified tables in `data/results/figure5/` |
| `results/figure6/` | Verified rewiring summary; full majority-edge tables in `data/results/figure6/` |
| `results/stimulation/` | Compact validation summaries; full response/cluster-audit tables in `data/results/stimulation/` |
| `docs/existing_raw_files.csv` | Original Zenodo raw-file URLs, sizes and checksums |
| Existing R Markdown/notebook files and `processed_data/` | Legacy materials retained for historical continuity; not certified as final-method reproductions |

Additional portable stimulation preparation/matching/statistics entrypoints are in [`code/stimulation/`](code/stimulation/README.md), with scope described in [reproduction coverage](docs/reproduction.md). To use saved inputs, decompress `data/upstream_inputs/*.csv.gz` to a separate project's `data/` and copy saved `*matches.csv` files into its `mario_out/`. Do not rerun matching merely to reconstruct the saved fused tables.

## Download and verify saved data

Large processed matrices, matching maps and fitted trajectories are distributed separately in the [published dataset v1.1.0](https://doi.org/10.5281/zenodo.23676187), not in GitHub release assets. Download **all six** derived-data archives:

```text
qbc_derived_data_v1.1.0_v6_part1of6.tar.gz
qbc_derived_data_v1.1.0_v6_part2of6.tar.gz
qbc_derived_data_v1.1.0_v6_part3of6.tar.gz
qbc_derived_data_v1.1.0_v6_part4of6.tar.gz
qbc_derived_data_v1.1.0_v6_part5of6.tar.gz
qbc_derived_data_v1.1.0_v6_part6of6.tar.gz
```

Place them in a `downloads/` directory under this repository root, together with the deposit's `README.txt`, `data_manifest.txt`, and `SHA256SUMS.txt`. Verify the eight outer checksum entries (six archives plus those two documents), then extract every part into this **same repository root**:

```bash
(cd downloads && sha256sum -c SHA256SUMS.txt)
for i in 1 2 3 4 5 6; do
  tar -xzf "downloads/qbc_derived_data_v1.1.0_v6_part${i}of6.tar.gz" -C .
done
(cd data && sha256sum -c SHA256SUMS)
```

These are normal independent archives, not binary fragments: do **not** concatenate them or create separate part-named extraction folders. Extraction order does not matter; all six are needed for the complete `data/`. Part 1 contains the global internal checksum list. Extract into a clean data location rather than mixing a different dataset version into an existing `data/`.

The repository's own `SHA256SUMS` is a separate curated code/documentation/table manifest. Zenodo's outer `SHA256SUMS.txt` does not hash the 19 retained original source archives; their historical MD5 provenance is documented separately in the deposit manifest and `docs/existing_raw_files.csv`.

Zenodo pins the exact analysis-source snapshot [`f0050c36f9f24ccad25b23c52b259c31346df92c`](https://github.com/aaron-tantk/QBC-CyTOF_cross-species_PBMCs/tree/f0050c36f9f24ccad25b23c52b259c31346df92c). The subsequent GitHub finalization changes documentation, citation metadata and related checksums, not scientific code or result bytes. The [original dataset DOI](https://doi.org/10.5281/zenodo.23406825) identifies version 1.0.0, not the v1.1.0 additions.

## Reproduce saved Figure 5 tables without refitting

Use Python with NumPy, pandas and SciPy. From the repository root:

```bash
python code/figure5/prepare_source.py --processed data/processed/primary --output work/source
python code/figure5/analysis_gmcsf.py prepare --source work/source --work work
```

Copy `data/figure5/graphs/` to `work/graphs/`, then:

```bash
python code/figure5/analysis_gmcsf.py summarize --work work
```

Compare `work/summary/` with `data/results/figure5/`. The preparation step excludes CD40L and recomputes the nine-marker proliferation score across Basal+GM-CSF rows within species. No new significance tests are performed. pSTAT5 is separate from the score.

For graph refitting and the distinction between reproducible summaries and stochastic refits, see [reproduction instructions](docs/reproduction.md) and [environment requirements](docs/environment.md).

## Reproduce saved Figure 6 summaries

```bash
python code/figure6/summarize_networks.py \
  --source data/figure6/bootstrap_support --output work/network_summary
```

Compare with `data/results/figure6/`. Majority directions summarize bootstrap output; they are not experimentally validated causal directions or a guarantee that the displayed edge set is acyclic.

## Licensing

Code: MIT, see `LICENSE`. Data and numerical result tables: CC-BY 4.0, see `DATA_LICENSE.md`. Original third-party datasets retain their attribution requirements. No manuscript, rendered figures or editorial revision notes are included in these additions.
