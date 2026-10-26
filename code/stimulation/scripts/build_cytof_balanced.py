"""Rebuild the balanced Bjornson-Hooper CyTOF dataset (MARIO input) from the
Zenodo CombinedMatrix parts.

- Reads 6 row-wise parts (732k cells, 183 donor x condition samples)
- Maps isotope channels -> marker names via Table S1 panel (mass-based)
- arcsinh(x/5) transformation (per Seurat-CyTOF-clustering.Rmd, cofactor=5)
- Balanced subsample per donor x condition: 260 cells (human) / 800 (rhesus),
  reproducing the original 35,880 / 36,000-cell balanced sets (seed 42)
- Writes cytof_{species}_balanced.csv with marker columns + donor + stimulation
"""
import re
import numpy as np
import pandas as pd

from portable import parser, paths, neutral_donors

PARTS = []
OUT_DIR = None

# mass -> marker name (fused-object naming convention, minus '-CyTOF' suffix)
MASS2MARKER = {
    113: ("CD235a", "CD233"), 115: ("CD45",)*2, 139: ("pSTAT3.Y705.",)*2, 140: ("CD61",)*2,
    141: ("CD7",)*2, 142: ("CD33",)*2, 143: ("CD11c",)*2,
    144: ("pMAPKAPK2.T334.",)*2, 145: ("pCREB.S133.",)*2, 146: ("pPLCg2.Y579.",)*2,
    147: ("pSTAT1.Y701.",)*2, 148: ("CD123",)*2, 149: ("pSTAT5.Y694.",)*2,
    150: ("pS6.S235.S236.",)*2, 151: ("CD14",)*2, 152: ("pERK1.2.T202.Y204.",)*2,
    153: ("CD11b",)*2, 154: ("pP38.MAPK.T180.Y182.",)*2, 155: ("CD8",)*2,
    156: ("CD4",)*2, 157: ("CD3",)*2, 158: ("CD66",)*2, 159: ("CD16",)*2,
    160: ("pZAP70.Syk.Y319.Y352.",)*2, 161: ("pTBK1.NAK.S172.",)*2,
    162: ("CD1c",)*2, 163: ("IkBa.amino.terminal",)*2, 164: ("BDCA3",)*2,
    165: ("pSTAT6.Y691.",)*2, 166: ("CD45RA",)*2, 167: ("FoxP3",)*2,
    168: ("CD161",)*2, 169: ("Ki67",)*2, 170: ("pSTAT4.Y693.",)*2,
    171: ("CCR7",)*2, 172: ("p4EBP1.T37.T46.",)*2, 173: ("CD19", "CD20"),
    174: ("IgM",)*2, 175: ("CD56",)*2, 176: ("HLA.DR",)*2,
}  # tuples: (human_name, rhesus_name)

def load_parts():
    dfs = []
    for i, p in enumerate(PARTS):
        df = pd.read_csv(p)
        dfs.append(df)
        print(f"part {i+1}: {df.shape}")
    full = pd.concat(dfs, ignore_index=True)
    print("combined:", full.shape)
    return full

def main():
    full = load_parts()
    print("Species counts:", full["Species"].value_counts().to_dict())
    print("Condition counts:", full["Condition"].value_counts().to_dict())

    # identify isotope columns and map to masses
    iso_cols = [c for c in full.columns if re.match(r"^[A-Z][a-z]?\d+Di$", c)]
    col_mass = {c: int(re.search(r"(\d+)", c).group(1)) for c in iso_cols}
    masses = sorted(col_mass.values())
    print(f"{len(iso_cols)} isotope cols; masses with panel mapping: {sorted(set(masses) & set(MASS2MARKER))}")
    missing = set(masses) - set(MASS2MARKER)
    print("unmapped masses (barcoding/other):", sorted(missing))

    frames = {}
    for sp_key, sp_name in [("Hu", "human"), ("Rh", "rhesus")]:
        sub = full[full["Species"] == sp_key].copy()
        idx = 0 if sp_key == "Hu" else 1
        rename = {c: MASS2MARKER[m][idx] for c, m in col_mass.items() if m in MASS2MARKER}
        markers = [rename[c] for c in sorted(rename, key=lambda c: col_mass[c])]
        sub = sub.rename(columns=rename)
        # arcsinh(x/5)
        sub[markers] = np.arcsinh(sub[markers] / 5.0)
        # balanced subsample per donor x condition
        cap = 260 if sp_key == "Hu" else 800
        rng = np.random.RandomState(42)
        keep = (sub.groupby(["Donor_Name", "Condition"], group_keys=False)
                  .apply(lambda g: g.sample(n=min(cap, len(g)), random_state=rng)))
        keep = keep.sample(frac=1.0, random_state=42).reset_index(drop=True)  # shuffle rows
        out = keep[markers + ["Donor_Name", "Condition"]].rename(
            columns={"Donor_Name": "donor", "Condition": "stimulation"})
        # Relabel only AFTER source group order and seeded selections are complete.
        out["donor"] = neutral_donors(out["donor"], sp_name, "cytof")
        out_path = f"{OUT_DIR}/cytof_{sp_name}_balanced.csv"
        out.to_csv(out_path, index=False)
        frames[sp_name] = out
        print(f"{sp_name}: {out.shape[0]} cells x {len(markers)} markers -> {out_path}")
        print("  per-condition:", out["stimulation"].value_counts().to_dict())
        print("  n donors:", out["donor"].nunique())
    return frames

if __name__ == "__main__":
    p = parser(__doc__, "raw/cytof", "data")
    args = paths(p, p.parse_args())
    PARTS = [args.input_dir / f"Combined_Data_part{i}of6.csv" for i in range(1, 7)]
    if any(not f.is_file() for f in PARTS):
        p.error("input directory must contain all six Combined_Data_partNof6.csv files")
    OUT_DIR = args.output_dir
    main()
