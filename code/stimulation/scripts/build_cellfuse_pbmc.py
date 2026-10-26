"""Build PBMC-only CellFuse reference CSVs (comment 3 design): exclude
Granulocyte/Myeloid, Platelets, Erythrocytes, Unresolved from the CyTOF reference.
Query (QBC) is unchanged - it is already PBMC."""
import os
import pandas as pd
from sklearn.model_selection import train_test_split

from portable import parser, paths, aligned_annotations, species_arg, species_list

p = parser(__doc__, "data", "cellfuse")
species_arg(p)
args = paths(p, p.parse_args())
DATA, OUT = args.input_dir, args.output_dir
os.makedirs(f"{OUT}/data", exist_ok=True)

ANCHORS = {
    "human": ["CD11b", "CD11c", "CD123", "CD14", "CD16", "CD161", "CD19", "CD1c",
              "CD3", "CD33", "CD4", "CD45", "CD45RA", "CD56", "CD7", "CD8",
              "HLA-DR", "IgM"],
    "rhesus": ["CD11b", "CD11c", "CD123", "CD14", "CD16", "CD161", "CD20", "CD1c",
               "CD3", "CD33", "CD4", "CD45", "CD45RA", "CD56", "CD7", "CD8",
               "HLA-DR", "IgM"],
}
EXCLUDE = {"Granulocyte/Myeloid", "Platelets", "Erythrocytes", "Unresolved"}

for sp in species_list(args):
    anchors = ANCHORS[sp]
    cyto = pd.read_csv(f"{DATA}/cytof_{sp}_balanced.csv")
    annot = pd.read_csv(f"{DATA}/cytof_{sp}_balanced_annotations.csv")
    aligned_annotations(cyto, annot, sp)
    assert (cyto["donor"].astype(str).values == annot["donor"].astype(str).values).all()
    cyto = cyto.rename(columns={"HLA.DR": "HLA-DR"})
    keep = ~annot["broad_lineage"].isin(EXCLUDE)
    ref = cyto.loc[keep.values, anchors].copy()
    ref["cluster.orig"] = annot.loc[keep.values, "broad_lineage"].values
    tr, val = train_test_split(ref, test_size=0.3, random_state=42,
                               stratify=ref["cluster.orig"])
    tr.to_csv(f"{OUT}/data/{sp}_cytof_pbmc_train.csv", index=False)
    val.to_csv(f"{OUT}/data/{sp}_cytof_pbmc_val.csv", index=False)
    print(f"{sp}: PBMC-only ref {ref.shape[0]} cells, "
          f"classes ({ref['cluster.orig'].nunique()}): "
          f"{ref['cluster.orig'].value_counts().to_dict()}")
