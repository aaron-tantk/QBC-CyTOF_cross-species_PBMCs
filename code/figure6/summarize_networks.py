"""Summarize archived bootstrap rows without fitting or changing support."""
from pathlib import Path
import argparse
import json
import pandas as pd

POPULATIONS = {
    "human": ["B_cells", "CD4_T_cells", "CD8_T_cells", "DC", "Monocytes", "NK_cells"],
    "rhesus": ["B_cells", "CD4_T_cells", "CD8_T_cells", "DC", "NK_cells", "T_activated"],
}
CONDITIONS = ["Basal", "CD40L", "GM-CSF"]


def consensus(path):
    df = pd.read_csv(path)
    assert set(["from", "to", "strength", "direction"]).issubset(df.columns)
    # Preserve archived rows: four pairs serialize at exactly 0.85 despite
    # the inference source's strict >0.85 comparison. Do not refilter them.
    assert (df.strength >= 0.85).all()
    assert df.direction.between(0, 1).all()
    df["pair"] = df.apply(lambda r: tuple(sorted([r["from"], r["to"]])), axis=1)
    rows = []
    for _, group in df.groupby("pair", sort=True):
        # Same first-maximum tie rule as the saved Figure 6 summary source.
        row = group.loc[group.direction.idxmax()]
        rows.append({k: row[k] for k in ["from", "to", "strength", "direction"]})
    return pd.DataFrame(rows, columns=["from", "to", "strength", "direction"])


def summarize(source, output):
    source, output = Path(source), Path(output)
    output.mkdir(parents=True, exist_ok=True)
    rewiring, inventory = [], []
    for species, populations in POPULATIONS.items():
        for population in populations:
            nets = {}
            for condition in CONDITIONS:
                name = f"fig6_arcs_{species}_{population}_{condition}.csv"
                result = consensus(source / name)
                result.to_csv(output / name.replace("fig6_arcs_", "majority_edges_"), index=False)
                nets[condition] = set(zip(result["from"], result["to"]))
                inventory.append(dict(species=species, population=population,
                                      condition=condition, n_majority_edges=len(result)))
            basal, gm = nets["Basal"], nets["GM-CSF"]
            # Compare unordered pairs so reversals count once, not gained+lost.
            b = {frozenset(edge): edge for edge in basal if any("pSTAT5" in x for x in edge)}
            g = {frozenset(edge): edge for edge in gm if any("pSTAT5" in x for x in edge)}
            gained, lost = len(g.keys()-b.keys()), len(b.keys()-g.keys())
            reversed_count = sum(b[k] != g[k] for k in b.keys() & g.keys())
            rewiring.append(dict(species=species, population=population,
                gained=gained, lost=lost, reversed=reversed_count,
                total=gained+lost+reversed_count, edges_basal=len(basal), edges_gmcsf=len(gm)))
    summary = pd.DataFrame(rewiring)
    summary.to_csv(output / "fig6_rewiring_summary.csv", index=False)
    pd.DataFrame(inventory).to_csv(output / "network_inventory.csv", index=False)
    print(summary.to_string(index=False))
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    summarize(args.source, args.output)
