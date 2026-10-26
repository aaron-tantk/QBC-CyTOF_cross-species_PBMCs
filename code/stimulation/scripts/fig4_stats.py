"""Re-derive Fig 4 intracellular-response statistics on saved fused datasets.

Design (source-label implementation described in manuscript v6):
- Per species x source broad population x intracellular node x stimulus:
  donor-mean delta (stim - basal) across CyTOF donors, paired t-test,
  requiring >= MIN_DONORS donors each with >= MIN_CELLS cells per condition.
- BH-FDR across eligible tests actually run (both species combined):
  375 primary or 480 whole-blood-reference sensitivity tests in saved results.
- Cross-species: Spearman rho of species mean donor deltas across nodes within
  each shared lineage, and globally across shared lineage x node pairs.
  Permutation p-values (10k): permute rhesus node deltas within lineage globally,
  and across nodes within a single lineage. Species labels are NOT swapped.

The 15 nodes are 14 phospho epitopes plus total IkBa.amino.terminal (IκBα).
Use --tag pbmc for the primary design; no tag selects whole-blood sensitivity.
Use --help for input/output paths. Source fine/broad per-cell labels are retained.
"""
import numpy as np
import pandas as pd
from scipy import stats
from itertools import product

from portable import parser, paths, neutral_donors

FUSED = RES = None

# Historical name retained: 14 phospho epitopes plus total IκBα.
PHOSPHO = ["pSTAT1.Y701.", "pSTAT3.Y705.", "IkBa.amino.terminal", "p4EBP1.T37.T46.",
           "pZAP70.Syk.Y319.Y352.", "pS6.S235.S236.", "pCREB.S133.", "pSTAT4.Y693.",
           "pSTAT5.Y694.", "pTBK1.NAK.S172.", "pSTAT6.Y691.", "pPLCg2.Y579.",
           "pERK1.2.T202.Y204.", "pP38.MAPK.T180.Y182.", "pMAPKAPK2.T334."]
STIMS = ["GM-CSF", "CD40L"]
MIN_DONORS = 5
MIN_CELLS = 3


def load_species(species, tag):
    suffix = f"_{tag}" if tag else ""
    f = pd.read_csv(f"{FUSED}/fused_{species}{suffix}.csv")
    phos_cols = [f"{p}-CyTOF" for p in PHOSPHO]
    assert all(c in f.columns for c in phos_cols), \
        f"{species}: missing {[c for c in phos_cols if c not in f.columns]}"
    df = f[["CyTOF_donor", "stimulation", "CyTOF_broad_lineage"] + phos_cols].copy()
    df.columns = ["donor", "stimulation", "population"] + PHOSPHO
    df["donor"] = neutral_donors(df["donor"], species, "cytof")
    return df


def donor_mean_deltas(df, population, node, stim):
    d = df[df["population"] == population]
    per = (d.groupby(["donor", "stimulation"])[node]
           .mean().unstack("stimulation"))
    if "Basal" not in per.columns or stim not in per.columns:
        return None
    n_basal = d[d.stimulation == "Basal"].groupby("donor")[node].size()
    n_stim = d[d.stimulation == stim].groupby("donor")[node].size()
    ok = per.index[(n_basal.reindex(per.index, fill_value=0) >= MIN_CELLS) &
                   (n_stim.reindex(per.index, fill_value=0) >= MIN_CELLS)]
    per = per.loc[ok].dropna(subset=["Basal", stim])
    if len(per) < MIN_DONORS:
        return None
    delta = per[stim] - per["Basal"]
    t, p = stats.ttest_rel(per[stim], per["Basal"])
    return {"n_donors": len(per), "mean_delta": delta.mean(),
            "median_delta": delta.median(), "t": t, "p": p,
            "deltas": delta}


def run(tag=""):
    suffix = f"_{tag}" if tag else ""
    data = {sp: load_species(sp, tag) for sp in ["human", "rhesus"]}
    pops = {sp: sorted(d["population"].unique()) for sp, d in data.items()}
    print({sp: (len(d), p) for sp, (d, p) in
           zip(["human", "rhesus"], [(data[s], pops[s]) for s in ["human", "rhesus"]])})

    # --- per-species response tests ---
    rows, delta_store = [], {}
    for sp in ["human", "rhesus"]:
        for pop in pops[sp]:
            for node in PHOSPHO:
                for stim in STIMS:
                    r = donor_mean_deltas(data[sp], pop, node, stim)
                    if r is None:
                        continue
                    rows.append({"species": sp, "population": pop, "node": node,
                                 "stim": stim, "n_donors": r["n_donors"],
                                 "mean_delta": r["mean_delta"],
                                 "median_delta": r["median_delta"], "p": r["p"]})
                    delta_store[(sp, pop, node, stim)] = r["deltas"]
    if not rows:
        raise ValueError("No eligible tests: require >=5 donors and >=3 cells per donor per condition")
    res = pd.DataFrame(rows)
    res["padj"] = stats.false_discovery_control(res["p"].values, method="bh") \
        if hasattr(stats, "false_discovery_control") else None
    if res["padj"].isnull().all():
        from statsmodels.stats.multitest import multipletests
        res["padj"] = multipletests(res["p"].values, method="fdr_bh")[1]
    res = res.sort_values("padj")
    res.to_csv(f"{RES}/fig4_phospho_stats{suffix}.csv", index=False)
    print(f"\ntests run: {len(res)} | FDR<0.05: {(res['padj'] < 0.05).sum()}")
    print("\ntop GM-CSF pSTAT5 responses:")
    sub = res[(res.node == "pSTAT5.Y694.") & (res.stim == "GM-CSF")].sort_values("mean_delta", ascending=False)
    print(sub[["species", "population", "n_donors", "mean_delta", "padj"]].to_string(index=False))

    # --- cross-species correlations ---
    shared = sorted(set(pops["human"]) & set(pops["rhesus"]))
    print(f"\nshared lineages: {shared}")
    xs_rows = []
    rng = np.random.default_rng(42)
    for stim in STIMS:
        # global: all shared-lineage x node delta pairs
        h, r = [], []
        for pop in shared:
            for node in PHOSPHO:
                k = (pop, node, stim)
                hh = delta_store.get(("human", pop, node, stim))
                rr = delta_store.get(("rhesus", pop, node, stim))
                if hh is not None and rr is not None:
                    h.append(hh.mean()); r.append(rr.mean())
        if len(h) >= 5:
            rho, p_asym = stats.spearmanr(h, r)
            # Stratified permutation: shuffle rhesus node deltas within each lineage.
            h_arr, r_arr = np.array(h), np.array(r)
            lin_arr = np.array([pop for pop in shared for node in PHOSPHO
                                if delta_store.get(("human", pop, node, stim)) is not None
                                and delta_store.get(("rhesus", pop, node, stim)) is not None])
            null = []
            for _ in range(10000):
                rn = r_arr.copy()
                for pop in np.unique(lin_arr):
                    idx = lin_arr == pop
                    rn[idx] = rng.permutation(rn[idx])
                null.append(stats.spearmanr(h_arr, rn).statistic)
            null = np.array(null)
            p_perm = (np.sum(np.abs(null) >= abs(rho)) + 1) / (len(null) + 1)
            xs_rows.append({"stim": stim, "scope": "global_shared_lineages",
                            "n_pairs": len(h), "rho": rho,
                            "p_perm": p_perm, "p_asym": p_asym})
            print(f"{stim} global cross-species rho={rho:.3f} (n={len(h)} lineage x node pairs, "
                  f"stratified perm p={p_perm:.4f}, asym p={p_asym:.2e})")
        # Per-lineage: rho across eligible intracellular nodes (up to 15).
        for pop in shared:
            h, r = [], []
            for node in PHOSPHO:
                hh = delta_store.get(("human", pop, node, stim))
                rr = delta_store.get(("rhesus", pop, node, stim))
                if hh is not None and rr is not None:
                    h.append(hh.mean()); r.append(rr.mean())
            if len(h) >= 8:
                rho, p_asym = stats.spearmanr(h, r)
                h_arr, r_arr = np.array(h), np.array(r)
                null = [stats.spearmanr(h_arr, rng.permutation(r_arr)).statistic
                        for _ in range(10000)]
                null = np.array(null)
                p_perm = (np.sum(np.abs(null) >= abs(rho)) + 1) / (len(null) + 1)
                xs_rows.append({"stim": stim, "scope": f"lineage:{pop}",
                                "n_pairs": len(h), "rho": rho,
                                "p_perm": p_perm, "p_asym": p_asym})
                print(f"  {stim} {pop}: rho={rho:.3f} across {len(h)} nodes "
                      f"(perm p={p_perm:.4f}, asym p={p_asym:.2e})")
    pd.DataFrame(xs_rows).to_csv(f"{RES}/fig4_cross_species{suffix}.csv", index=False)
    return res


if __name__ == "__main__":
    p = parser(__doc__, "fused", "results")
    p.add_argument("--tag", choices=("", "pbmc"), default="")
    args = paths(p, p.parse_args())
    FUSED, RES = args.input_dir, args.output_dir
    run(args.tag)
