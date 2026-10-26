"""Re-run MARIO matching (QBC 30k vs CyTOF balanced 30k) with the exact
parameters from mario-pY-CyTOF_stim-QBC-matching_v2.ipynb.

Use --help for explicit input/output paths. Matching is stochastic; this port
does not set a new global seed. --smoke is an OPTIONAL stochastic rerun, not a
release verification step.
"""
import os, sys, json, time
import numpy as np
import pandas as pd

from portable import parser, paths, neutral_donors, species_arg, species_list

DATA = None
OUT = None

ANCHORS_HU = ["CD11b","CD11c","CD123","CD14","CD16","CD161","CD19","CD1c","CD3","CD33",
              "CD4","CD45","CD45RA","CD56","CD7","CD8","HLA-DR","IgM"]
ANCHORS_RH = ["CD11b","CD11c","CD123","CD14","CD16","CD161","CD20","CD1c","CD3","CD33",
              "CD4","CD45","CD45RA","CD56","CD7","CD8","HLA-DR","IgM"]

def load_species(species, smoke=False, pbmc_only=False):
    qbc = pd.read_csv(f"{DATA}/qbc_{species}_30k.csv")
    cyto_file = f"{DATA}/cytof_{species}_balanced{'_pbmc_only' if pbmc_only else ''}.csv"
    cyto = pd.read_csv(cyto_file)
    qbc["donor"] = neutral_donors(qbc["donor"], species, "qbc")
    cyto["donor"] = neutral_donors(cyto["donor"], species, "cytof")
    anchors = ANCHORS_HU if species == "human" else ANCHORS_RH
    # MARIO convention: overlapping features must have identical names in both datasets.
    # Rename CyTOF cols to match QBC where they are anchors; keep others distinct.
    cyto = cyto.rename(columns={"HLA.DR": "HLA-DR"})
    # shuffle rows (MARIO best practice; QBC csv is donor-ordered)
    qbc = qbc.sample(frac=1.0, random_state=42).reset_index(drop=True)
    cyto = cyto.sample(frac=1.0, random_state=42).reset_index(drop=True)
    if smoke:
        qbc = qbc.iloc[:3000].reset_index(drop=True)
        cyto = cyto.iloc[:3000].reset_index(drop=True)
    elif not pbmc_only:
        # notebook subsampled CyTOF 35.9k/36k -> 30k (seed 42); QBC already 30k
        cyto = cyto.iloc[:30000].reset_index(drop=True)
    return qbc, cyto, anchors

def run_species(species, smoke=False, pbmc_only=False):
    from mario.match import pipelined_mario
    t0 = time.time()
    qbc, cyto, anchors = load_species(species, smoke, pbmc_only)
    print(f"[{species}] QBC {qbc.shape}, CyTOF {cyto.shape}; anchors: {len(anchors)}")

    # MARIO-ready arrays: all features except metadata; shared anchors aligned by name.
    # pipelined_mario takes dataframes with identical column order for shared features?
    # Notebook passed numpy arrays with ALL features (134 and 40) and MARIO aligns by
    # position of overlapping features -- it requires shared features to be identified
    # by name, so pass DataFrames (mario-py accepts DataFrames and uses column names).
    qbc_X = qbc.drop(columns=["donor"])
    cyto_X = cyto.drop(columns=["donor", "stimulation"])
    # keep only shared anchors as overlapping; make non-shared names unique per modality
    qbc_X = qbc_X.rename(columns={c: c for c in qbc_X.columns})  # QBC names as-is
    cyto_X = cyto_X.rename(columns={c: (c if c in anchors else c + "_CY") for c in cyto_X.columns})

    kw = dict(normalization=True, n_batches=5, n_matched_per_cell=1,
              sparsity_ovlp=400, sparsity_all=400,
              n_components_ovlp=10, n_components_all=12,
              n_cancor=5, n_wts=4, n_clusters=10, n_components_filter=10,
              bad_prop=0.2, max_iter_filter=20, knn=False,
              embed_dim=10, max_iter_embed=(50 if smoke else 500),
              save_path=OUT, verbose=True)
    matching_lst, embedding_lst = pipelined_mario(data_lst=[qbc_X, cyto_X], **kw)

    # pipelined_mario writes save_path/final_matching.csv (row_idx_0, row_idx_1; NA = unmatched)
    tag = ("smoke_" if smoke else "") + ("pbmc_" if pbmc_only else "")
    m = pd.read_csv(f"{OUT}/final_matching.csv")
    os.rename(f"{OUT}/final_matching.csv", f"{OUT}/{tag}{species}_final_matching_raw.csv")
    for i in range(2):
        if os.path.exists(f"{OUT}/embedding_{i}.csv"):
            os.rename(f"{OUT}/embedding_{i}.csv", f"{OUT}/{tag}{species}_embedding_{i}.csv")
    m = m.dropna()
    m["row_idx_1"] = m["row_idx_1"].astype(int)
    matches = m[["row_idx_0", "row_idx_1"]].values  # (n_pairs, 2): QBC idx, CyTOF idx
    np.save(f"{OUT}/{tag}{species}_matches.npy", matches)
    m.to_csv(f"{OUT}/{tag}{species}_matches.csv", index=False)
    meta = {"species": species, "n_pairs": int(matches.shape[0]),
            "runtime_min": round((time.time()-t0)/60, 2), "smoke": smoke,
            "params": {k: (v if isinstance(v, (int, float, bool, str)) else str(v)) for k, v in kw.items()}}
    with open(f"{OUT}/{tag}{species}_run_meta.json", "w") as f:
        json.dump(meta, f, indent=2)
    print(f"[{species}] matched pairs: {matches.shape[0]} in {meta['runtime_min']} min")
    return matches, qbc, cyto

if __name__ == "__main__":
    p = parser(__doc__, "data", "mario_out")
    species_arg(p)
    p.add_argument("--smoke", action="store_true")
    p.add_argument("--pbmc-only", action="store_true")
    p.add_argument("--prepare-only", action="store_true",
                   help="Export exact shuffled MARIO inputs with saved source labels; do not match")
    args = paths(p, p.parse_args())
    DATA, OUT = args.input_dir, str(args.output_dir)
    for sp in species_list(args):
        if args.prepare_only:
            from portable import aligned_annotations
            bal = pd.read_csv(f"{DATA}/cytof_{sp}_balanced.csv")
            ann = pd.read_csv(f"{DATA}/cytof_{sp}_balanced_annotations.csv")
            aligned_annotations(bal, ann, sp)
            if args.pbmc_only:
                keep = ~ann["broad_lineage"].isin(
                    {"Granulocyte/Myeloid", "Platelets", "Erythrocytes", "Unresolved"})
                selected = pd.read_csv(f"{DATA}/cytof_{sp}_balanced_pbmc_only.csv")
                feats = [c for c in bal if c not in ("donor", "stimulation")]
                if len(selected) != keep.sum() or not np.allclose(
                        selected[feats].to_numpy(), bal.loc[keep, feats].to_numpy()):
                    raise ValueError(f"{sp}: PBMC input does not match source label filter")
                ann = ann.loc[keep].reset_index(drop=True)
            q, c, _ = load_species(sp, args.smoke, args.pbmc_only)
            ann = ann.sample(frac=1.0, random_state=42).reset_index(drop=True)
            ann = ann.iloc[:len(c)]
            c["CyTOF_annotation"] = ann["CyTOF_annotation"].to_numpy()
            c["broad_lineage"] = ann["broad_lineage"].to_numpy()
            prefix = ("smoke_" if args.smoke else "") + ("pbmc_" if args.pbmc_only else "")
            q.to_csv(f"{OUT}/{prefix}mario_input_qbc_{sp}.csv", index=False)
            c.to_csv(f"{OUT}/{prefix}mario_input_cytof_{sp}.csv", index=False)
        else:
            run_species(sp, args.smoke, args.pbmc_only)
