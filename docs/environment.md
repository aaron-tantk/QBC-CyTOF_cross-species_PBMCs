# Environment and reproducibility boundary

## Table reconstruction

Python: NumPy, pandas, SciPy. Optional plotting source also requires matplotlib and adjustText; artwork is not included in these additions.

Exact installed versions used for the earlier source-package validation are recorded in `validation_summary.json`. This is a historical observed validation environment, not a reconstructed lockfile for every historical analysis or a new validation run during documentation finalization.

## Figure 5 graph refitting

Required: Monocle3 **1.3.1**, igraph **1.6.0**, data.table, jsonlite, RhpcBLASctl and their R dependencies. The saved fitter enforces the first two versions and requires a Unix execution environment. Original seed handling and defaults are preserved.

## Other recovered sources

Seurat, flowCore (legacy FCS import), bnlearn (optional network refit), Harmony and MARIO dependencies are needed only for the corresponding optional stages. Refer to upstream instructions for recovered implementation details. Do not treat currently installed versions as proven historical versions. No new MARIO, clustering or network bootstrap was run for source packaging or documentation finalization.

Generated summaries can be compared against the saved numerical tables without running these stochastic stages. A new environment can change refitted embeddings, clusters and network structures even when parameters are nominally identical.
