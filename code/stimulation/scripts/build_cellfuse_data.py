"""Build CellFuse-format CSVs for QBC->CyTOF cross-modality label transfer.

Reference = annotated balanced CyTOF (18 shared anchors + broad_lineage label),
split 70/30 stratified (per CellFuse DataPreparation guide).
Query = QBC 30k cells, same 18 anchors in identical column order.
"""
import os
import numpy as np
import pandas as pd
from sklearn.model_selection import train_test_split

from portable import parser, paths, aligned_annotations, species_arg, species_list

p = parser(__doc__, "data", "cellfuse")
species_arg(p)
args = paths(p, p.parse_args())
DATA, OUT = args.input_dir, args.output_dir
os.makedirs(f"{OUT}/data", exist_ok=True)
os.makedirs(f"{OUT}/Saved_model", exist_ok=True)
os.makedirs(f"{OUT}/pred", exist_ok=True)

ANCHORS = {
    "human": ["CD11b", "CD11c", "CD123", "CD14", "CD16", "CD161", "CD19", "CD1c",
              "CD3", "CD33", "CD4", "CD45", "CD45RA", "CD56", "CD7", "CD8",
              "HLA-DR", "IgM"],
    "rhesus": ["CD11b", "CD11c", "CD123", "CD14", "CD16", "CD161", "CD20", "CD1c",
               "CD3", "CD33", "CD4", "CD45", "CD45RA", "CD56", "CD7", "CD8",
               "HLA-DR", "IgM"],
}

for sp in species_list(args):
    anchors = ANCHORS[sp]
    # --- reference: CyTOF balanced + annotations ---
    cyto = pd.read_csv(f"{DATA}/cytof_{sp}_balanced.csv")
    annot = pd.read_csv(f"{DATA}/cytof_{sp}_balanced_annotations.csv")
    aligned_annotations(cyto, annot, sp)
    assert len(cyto) == len(annot), f"{sp}: row mismatch"
    # row-order alignment: donor+stimulation must match row-by-row
    assert (cyto["donor"].astype(str).values == annot["donor"].astype(str).values).all(), \
        f"{sp}: donor misalignment"
    assert (cyto["stimulation"].astype(str).values == annot["stimulation"].astype(str).values).all(), \
        f"{sp}: stimulation misalignment"
    cyto = cyto.rename(columns={"HLA.DR": "HLA-DR"})
    missing = [a for a in anchors if a not in cyto.columns]
    assert not missing, f"{sp} CyTOF missing {missing}"
    ref = cyto[anchors].copy()
    ref["cluster.orig"] = annot["broad_lineage"].values

    tr, val = train_test_split(ref, test_size=0.3, random_state=42,
                               stratify=ref["cluster.orig"])
    tr.to_csv(f"{OUT}/data/{sp}_cytof_train.csv", index=False)
    val.to_csv(f"{OUT}/data/{sp}_cytof_val.csv", index=False)

    # --- query: QBC 30k ---
    qbc = pd.read_csv(f"{DATA}/qbc_{sp}_30k.csv")
    qbc = qbc.rename(columns={"HLA.DR": "HLA-DR"})
    missing = [a for a in anchors if a not in qbc.columns]
    assert not missing, f"{sp} QBC missing {missing}"
    qbc[anchors].to_csv(f"{OUT}/data/{sp}_qbc_test.csv", index=False)

    print(f"{sp}: train {tr.shape}, val {val.shape}, test {qbc[anchors].shape}")
    print(f"  classes ({ref['cluster.orig'].nunique()}):",
          sorted(ref["cluster.orig"].unique()))
    print("  class counts:", ref["cluster.orig"].value_counts().to_dict())
