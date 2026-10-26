"""Assemble corrected fused QBC<->CyTOF paired datasets from the MARIO re-run.

Reproduces the exact row shuffle/head used in run_mario.py (deterministic, seed 42),
joins match indices back to unshuffled rows via row_id, attaches CyTOF annotations,
and writes fused_{species}.csv plus matching-QC reports.
"""
import numpy as np
import pandas as pd

from pathlib import Path
import os
from portable import parser, paths, neutral_donors, aligned_annotations, species_arg, species_list

DATA = OUT = FUSED = None

ANCHORS = {
    "human": ["CD11b","CD11c","CD123","CD14","CD16","CD161","CD19","CD1c","CD3","CD33",
              "CD4","CD45","CD45RA","CD56","CD7","CD8","HLA-DR","IgM"],
    "rhesus": ["CD11b","CD11c","CD123","CD14","CD16","CD161","CD20","CD1c","CD3","CD33",
               "CD4","CD45","CD45RA","CD56","CD7","CD8","HLA-DR","IgM"],
}

EXCLUDE = {"Granulocyte/Myeloid", "Platelets", "Erythrocytes", "Unresolved"}

def assemble(species, tag=""):
    pbmc = tag == "pbmc"
    # --- reload inputs exactly as run_mario.py did ---
    qbc = pd.read_csv(f"{DATA}/qbc_{species}_30k.csv")
    qbc["donor"] = neutral_donors(qbc["donor"], species, "qbc")
    ann = pd.read_csv(f"{DATA}/cytof_{species}_balanced_annotations.csv")
    bal = pd.read_csv(f"{DATA}/cytof_{species}_balanced.csv")
    aligned_annotations(bal, ann, species)
    if pbmc:
        cyto = pd.read_csv(f"{DATA}/cytof_{species}_balanced_pbmc_only.csv")
        # map pbmc_only rows back to balanced rows for annotation:
        # rebuild the keep-mask from balanced + annotations and verify identity
        keep_idx = np.flatnonzero(~ann["broad_lineage"].isin(EXCLUDE).values)
        marker_cols = [c for c in bal.columns if c not in ("donor", "stimulation")]
        assert len(keep_idx) == len(cyto), f"{species}: pbmc row count mismatch"
        assert np.allclose(bal[marker_cols].values[keep_idx], cyto[marker_cols].values), \
            f"{species}: pbmc_only rows do not match balanced[keep_idx]"
        cyto["donor"] = bal.iloc[keep_idx]["donor"].to_numpy()
        if not (cyto["stimulation"].astype(str).to_numpy() ==
                bal.iloc[keep_idx]["stimulation"].astype(str).to_numpy()).all():
            raise ValueError(f"{species}: PBMC stimulation order mismatch")
        cyto["row_id"] = keep_idx  # row_id = index into balanced/annotation tables
    else:
        cyto = bal.copy()
        cyto["row_id"] = np.arange(len(cyto))
    qbc["row_id"] = np.arange(len(qbc))
    qbc_s = qbc.sample(frac=1.0, random_state=42).reset_index(drop=True)
    cyto_s = cyto.sample(frac=1.0, random_state=42).reset_index(drop=True)
    if not pbmc:
        cyto_s = cyto_s.iloc[:30000].reset_index(drop=True)

    m = pd.read_csv(f"{OUT}/{tag + '_' if tag else ''}{species}_matches.csv")  # row_idx_0 -> qbc_s, row_idx_1 -> cyto_s
    for col, n in (("row_idx_0", len(qbc_s)), ("row_idx_1", len(cyto_s))):
        v = m[col]
        if v.isna().any() or (v != np.floor(v)).any() or ((v < 0) | (v >= n)).any():
            raise ValueError(f"{species}: invalid zero-based match index in {col}")
        m[col] = v.astype(int)
    q = qbc_s.loc[m["row_idx_0"]].reset_index(drop=True)
    c = cyto_s.loc[m["row_idx_1"]].reset_index(drop=True)

    # --- attach CyTOF annotations (row_id = unshuffled balanced-CSV row = annotation CSV row) ---
    c["CyTOF_annotation"] = ann.loc[c["row_id"], "CyTOF_annotation"].values
    c["broad_lineage"] = ann.loc[c["row_id"], "broad_lineage"].values

    # --- fused dataframe ---
    qbc_feats = [col for col in q.columns if col not in ("donor", "row_id")]
    cyto_feats = [col for col in c.columns if col not in ("donor", "stimulation", "row_id", "CyTOF_annotation", "broad_lineage")]
    fused = pd.concat([
        q[qbc_feats].add_suffix("-QBC").reset_index(drop=True),
        c[cyto_feats].add_suffix("-CyTOF").reset_index(drop=True),
        pd.DataFrame({"QBC_donor": q["donor"].values, "CyTOF_donor": c["donor"].values,
                      "stimulation": c["stimulation"].values,
                      "CyTOF_annotation": c["CyTOF_annotation"].values,
                      "CyTOF_broad_lineage": c["broad_lineage"].values,
                      "qbc_row_id": q["row_id"].values, "cytof_row_id": c["row_id"].values}),
    ], axis=1)
    # --- QC 1: within-pair shared-anchor correlation ---
    anchors = ANCHORS[species]
    Q = q[anchors].values
    # cyto csv keeps HLA.DR column name; map anchor HLA-DR -> HLA.DR
    C = c[[("HLA.DR" if a == "HLA-DR" else a) for a in anchors]].values
    r = np.array([np.corrcoef(Q[i], C[i])[0, 1] for i in range(len(Q))])
    fused["within_pair_anchor_r"] = r
    suffix = f"_{tag}" if tag else ""
    fused.to_csv(f"{FUSED}/fused_{species}{suffix}.csv", index=False)

    # --- QC 2: per-population match rates (CyTOF side, within the 30k MARIO input) ---
    cyto_s_ann = cyto_s.copy()
    cyto_s_ann["CyTOF_annotation"] = ann.loc[cyto_s_ann["row_id"], "CyTOF_annotation"].values
    cyto_s_ann["broad_lineage"] = ann.loc[cyto_s_ann["row_id"], "broad_lineage"].values
    cyto_s_ann["matched"] = False
    cyto_s_ann.loc[np.unique(m["row_idx_1"]), "matched"] = True
    rates = (cyto_s_ann.groupby("broad_lineage")["matched"]
             .agg(n_input="count", n_matched="sum", match_rate="mean")
             .sort_values("match_rate"))
    rates.to_csv(f"{FUSED}/{species}{suffix}_cytof_match_rates.csv")

    # QBC-side match rate per donor
    qbc_s2 = qbc_s.copy()
    qbc_s2["matched"] = False
    qbc_s2.loc[np.unique(m["row_idx_0"]), "matched"] = True
    qbc_rates = qbc_s2.groupby("donor")["matched"].agg(n_input="count", n_matched="sum", match_rate="mean")
    qbc_rates.to_csv(f"{FUSED}/{species}{suffix}_qbc_donor_match_rates.csv")

    print(f"\n=== {species} ===")
    print(f"pairs: {len(fused)} | within-pair anchor r: median={np.median(r):.3f}, "
          f"frac>0.5={100*(r>0.5).mean():.1f}%, frac>0.7={100*(r>0.7).mean():.1f}%")
    print("CyTOF per-lineage match rates:")
    print(rates.round(3).to_string())
    print("QBC per-donor match rates:")
    print(qbc_rates.round(3).to_string())
    print("stimulation distribution of pairs:", fused["stimulation"].value_counts().to_dict())
    return fused

if __name__ == "__main__":
    p = parser(__doc__, "data", "fused")
    species_arg(p)
    p.add_argument("--matches-dir", help="Overrides PROJECT_ROOT/mario_out")
    p.add_argument("--tag", choices=("", "pbmc"), default="")
    args = paths(p, p.parse_args())
    if not args.matches_dir and not args.project_root:
        p.error("--matches-dir or --project-root/QBC_PROJECT_ROOT is required")
    DATA, FUSED = args.input_dir, args.output_dir
    OUT = Path(args.matches_dir) if args.matches_dir else Path(args.project_root) / "mario_out"
    if not OUT.is_dir():
        p.error(f"matches directory does not exist: {OUT}")
    for sp in species_list(args):
        assemble(sp, args.tag)
