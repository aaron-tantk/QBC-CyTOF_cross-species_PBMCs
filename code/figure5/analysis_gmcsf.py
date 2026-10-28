"""GM-CSF-only Figure 5: exact row filtering and descriptive donor summaries."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

SPECIES = ["human", "rhesus"]
GROUPS = ["Basal", "GM-CSF"]
SIGNATURE = ["Ki67-CyTOF", "HLA.DR-CyTOF", "pS6.S235.S236.-CyTOF",
             "p4EBP1.T37.T46.-CyTOF", "pERK1.2.T202.Y204.-CyTOF",
             "HLA-DR-QBC", "CD71-QBC", "CD25-QBC", "CD69-QBC"]
STAT5 = "pSTAT5.Y694.-CyTOF"
TYPES = ["CD4_T_cells", "CD8_T_cells", "Activated_T_cells"]
COUNTS = {"human": {"Basal": 7740, "GM-CSF": 6958},
          "rhesus": {"Basal": 4807, "GM-CSF": 4989}}


def prepare(source, work):
    source, work = Path(source), Path(work)
    out = work / "inputs"
    out.mkdir(parents=True, exist_ok=True)
    audit = {}
    for sp in SPECIES:
        m = pd.read_csv(source / f"{sp}_metadata.csv.gz")
        x = pd.read_csv(source / f"{sp}_matrix.csv.gz", float_precision="round_trip")
        assert m.cell_id.equals(x.cell_id) and m.cell_id.is_unique
        assert set(m.condition) == {"Basal", "CD40L", "GM-CSF"}
        use = m.condition.isin(GROUPS)
        m, x = m.loc[use].reset_index(drop=True), x.loc[use].reset_index(drop=True)
        assert m.condition.value_counts().to_dict() == COUNTS[sp]
        assert m.cell_id.equals(x.cell_id)
        assert x.shape[1] == 175 and len([c for c in x if c.endswith("-QBC")]) == 134
        vals = x[SIGNATURE].to_numpy()
        # Representative pilot and full-data formula checks use identical columns.
        assert np.isfinite(vals[:100]).all() and np.isfinite(x.iloc[:, 1:].to_numpy()).all()
        mu, sd = vals.mean(axis=0), vals.std(axis=0, ddof=1)
        assert (sd > 0).all() and len(SIGNATURE) == 9 and STAT5 not in SIGNATURE
        score = ((vals - mu) / sd).mean(axis=1)
        np.testing.assert_allclose(score[:100], ((x[SIGNATURE].iloc[:100] - mu) / sd).mean(axis=1),
                                   rtol=1e-12, atol=1e-12)
        m["nine_marker_score"] = score
        m["reference_group"] = m.condition
        m["is_T"] = m.population.isin(TYPES)
        assert not m.condition.eq("CD40L").any()
        m.to_csv(out / f"{sp}_metadata.csv.gz", index=False)
        x.to_csv(out / f"{sp}_matrix.csv.gz", index=False)
        pd.DataFrame({"species": sp, "feature": SIGNATURE, "mean": mu, "sd": sd,
                      "scale_population": "Basal_and_GM-CSF_fused_rows_within_species"
                      }).to_csv(out / f"{sp}_signature_scale.csv", index=False)
        audit[sp] = {"condition_counts": COUNTS[sp], "n_features": 174,
                     "n_unique_references": int(m.reference_id.nunique()),
                     "max_reference_reuse": int(m.reference_id.value_counts().max()),
                     "signature": SIGNATURE, "CD40L_rows": 0,
                     "source_matrix_sha256": hashlib.sha256(
                         (source / f"{sp}_matrix.csv.gz").read_bytes()).hexdigest()}
    (out / "input_validation.json").write_text(json.dumps(audit, indent=2))
    print(json.dumps(audit, indent=2))
    return audit


def donor_summaries(m):
    reference = m.groupby(["species", "reference_id", "cytof_donor", "population",
                           "reference_group"], as_index=False).agg(
        nine_marker_score=("nine_marker_score", "mean"),
        pSTAT5=("pSTAT5", "mean"), n_links=("cell_id", "size"))
    donor = reference.groupby(["species", "cytof_donor", "population", "reference_group"],
                              as_index=False).agg(
        nine_marker_score=("nine_marker_score", "mean"), pSTAT5=("pSTAT5", "mean"),
        n_unique_reference_cells=("reference_id", "size"))
    changes, paired_all = [], []
    for (sp, pop), d in donor.groupby(["species", "population"]):
        counts = d.pivot(index="cytof_donor", columns="reference_group",
                         values="n_unique_reference_cells").reindex(columns=GROUPS).fillna(0)
        for measure in ["nine_marker_score", "pSTAT5"]:
            v = d.pivot(index="cytof_donor", columns="reference_group",
                        values=measure).reindex(columns=GROUPS)
            use = v.notna().all(axis=1) & counts.ge(3).all(axis=1)
            pairs = v.loc[use].copy()
            pairs["delta"] = pairs["GM-CSF"] - pairs["Basal"]
            pairs["n_basal_reference_cells"] = counts.loc[use, "Basal"].astype(int)
            pairs["n_GMCSF_reference_cells"] = counts.loc[use, "GM-CSF"].astype(int)
            pairs["species"], pairs["population"], pairs["measure"] = sp, pop, measure
            pairs["cytof_donor"] = pairs.index
            paired_all.append(pairs.reset_index(drop=True))
            z = pairs.delta
            changes.append({"species": sp, "population": pop, "measure": measure,
                            "n_paired_donors": len(z), "eligible_for_ranking": len(z) >= 5,
                            "mean_delta": z.mean(), "median_delta": z.median(),
                            "donor_sd": z.std(ddof=1), "donor_min": z.min(),
                            "donor_max": z.max()})
    change = pd.DataFrame(changes)
    group = ["species", "measure"]
    valid = change.eligible_for_ranking
    change.loc[valid, "rank_absolute_change"] = change.loc[valid].assign(
        effect=change.loc[valid, "mean_delta"].abs()).groupby(group).effect.rank(
            ascending=False, method="min")
    change.loc[valid, "rank_signed_change"] = change.loc[valid].groupby(group).mean_delta.rank(
        ascending=False, method="min")
    return donor, pd.concat(paired_all, ignore_index=True), change


def summarize(work):
    work = Path(work)
    out = work / "summary"
    out.mkdir(exist_ok=True)
    metadata, matrices, audit = [], {}, {}
    for sp in SPECIES:
        m = pd.read_csv(work / "inputs" / f"{sp}_metadata.csv.gz")
        x = pd.read_csv(work / "inputs" / f"{sp}_matrix.csv.gz", float_precision="round_trip")
        assert m.cell_id.equals(x.cell_id)
        assert m.condition.value_counts().to_dict() == COUNTS[sp]
        scale = pd.read_csv(work / "inputs" / f"{sp}_signature_scale.csv")
        independently_scaled = (x[SIGNATURE] - x[SIGNATURE].mean()) / x[SIGNATURE].std(ddof=1)
        np.testing.assert_allclose(independently_scaled.mean(axis=1), m.nine_marker_score,
                                   rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(scale["mean"], x[SIGNATURE].mean(), rtol=1e-12, atol=1e-12)
        for key in ["cytof_donor", "population", "condition"]:
            assert m.groupby("reference_id")[key].nunique().eq(1).all()
        m["pSTAT5"] = x[STAT5]
        metadata.append(m)
        matrices[sp] = x.set_index("cell_id")
    donor, paired, change = donor_summaries(pd.concat(metadata, ignore_index=True))
    donor.to_csv(out / "donor_population_scores_pSTAT5.csv", index=False)
    paired.to_csv(out / "paired_donor_changes.csv", index=False)
    change.to_csv(out / "population_changes.csv", index=False)
    allmarkers, coverage = [], []
    for sp in SPECIES:
        for cond in GROUPS:
            cells = pd.read_csv(work / "graphs/tables" / f"{sp}_{cond}_cells.csv.gz")
            assert cells.condition.eq(cond).all() and len(cells) == COUNTS[sp][cond]
            roots = cells.nine_marker_score >= cells.nine_marker_score.quantile(.95)
            np.testing.assert_array_equal(roots, cells.root_cell)
            np.testing.assert_array_equal(np.isfinite(cells.pseudotime), cells.finite_pseudotime)
            for pop, d in cells.groupby("population"):
                coverage.append({"species": sp, "condition": cond, "population": pop,
                                 "n_cells": len(d), "n_finite": int(d.finite_pseudotime.sum())})
        t = cells.loc[cells.is_T & cells.finite_pseudotime]
        assert len(t) > 100 and t.population.isin(TYPES).all()
        exp = matrices[sp].loc[t.cell_id]
        qcols = [c for c in exp if c.endswith("-QBC")]
        assert len(qcols) == 134
        for feature in qcols:
            rho = spearmanr(exp[feature], t.pseudotime).statistic
            allmarkers.append({"species": sp, "condition": "GM-CSF", "feature": feature,
                               "spearman_rho": rho, "abs_rho": abs(rho), "n_T_cells": len(t),
                               "included_in_signature": feature in SIGNATURE})
        audit[sp] = {"n_T_GMCSF": int(cells.is_T.sum()), "n_T_finite": len(t)}
    markers = pd.DataFrame(allmarkers).sort_values(
        ["species", "abs_rho", "feature"], ascending=[True, False, True])
    markers["rank_within_species"] = markers.groupby("species").cumcount() + 1
    markers.to_csv(out / "T_cell_QBC_marker_trends.csv", index=False)
    pd.DataFrame(coverage).to_csv(out / "population_pseudotime_coverage.csv", index=False)
    graphs = [json.loads(p.read_text()) for p in sorted(
        (work / "graphs/tables").glob("*_summary.json"))]
    assert len(graphs) == 4
    pd.DataFrame(graphs).to_csv(out / "graph_coverage.csv", index=False)
    summary = {"condition": "GM-CSF_vs_Basal_only", "T_cell_coverage": audit, "graphs": graphs,
               "top_markers": {sp: markers[markers.species == sp].head(20).to_dict("records")
                               for sp in SPECIES},
               "population_changes": change.replace({np.nan: None}).to_dict("records"),
               "significance_tests": False}
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    print(change.sort_values(["measure", "species", "rank_absolute_change"]).to_string(index=False))
    print("\nT-CELL COVERAGE", audit)
    print(markers.groupby("species").head(5).to_string(index=False))
    return summary


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("action", choices=["prepare", "summarize"])
    ap.add_argument("--work", type=Path, required=True)
    ap.add_argument("--source", type=Path)
    a = ap.parse_args()
    if a.action == "prepare":
        if a.source is None:
            ap.error("--source is required for prepare")
        prepare(a.source, a.work)
    else:
        summarize(a.work)
