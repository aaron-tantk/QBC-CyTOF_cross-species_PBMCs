# QBC-CyTOF cross-species PBMCs

Analysis code and processed marker tables for the cross-species (human & rhesus macaque)
PBMC study combining **Quantum Barcoding (QBC)** protein profiling with **CyTOF** mass
cytometry under basal and stimulated (GM-CSF, CD40L) conditions.

This repository accompanies the manuscript (Cell Systems, under review; preprint on bioRxiv).

## Repository contents

```
code/       Analysis notebooks / R Markdown scripts
markers/    Per-cluster marker tables and cluster-proportion statistics (CSV)
```

### `code/`
| File | Purpose |
|---|---|
| `Seurat-QBC-clustering.Rmd` | QBC surface-marker clustering (Seurat) |
| `Seurat-CyTOF-clustering.Rmd` | CyTOF clustering (Seurat) |
| `Seurat-CyTOF-QBC_fused_dataset_clustering.Rmd` | Fused CyTOF–QBC clustering |
| `Seurat-CyTOF-QBC_fused_dataset_clustering_stimulations.Rmd` | Fused clustering across stimulation conditions |
| `mario-pY-CyTOF_stim-QBC-matching_v2.ipynb` | MARIO integration / phospho-CyTOF–QBC matching |

### `markers/`
Top-20 marker tables per cluster (CyTOF-only and fused) and cluster-proportion
statistical tests for human and rhesus macaque.

## Processed objects & raw data

Processed Seurat objects (`.rds`) and the raw CyTOF/QBC data (~1.4 GB) are archived on
**Zenodo**: DOI *(to be added on publication)*.

## Requirements

- **R** (≥4.x): `Seurat`, `flowCore`, `monocle3`, `bnlearn`, `igraph`, `tidygraph`,
  `ggraph`, `ComplexHeatmap`, `pheatmap`, `ggplot2`, `RColorBrewer`
- **Python** (≥3.x) for the MARIO notebook: `mario`, `scanpy`/anndata, `numpy`, `pandas`

## License

Code in this repository is released under the **MIT License** (see `LICENSE`).
Figures and tables are distributed with the manuscript preprint (bioRxiv), not here.

## Citation

*(Manuscript citation and Zenodo DOI to be added on publication.)*
