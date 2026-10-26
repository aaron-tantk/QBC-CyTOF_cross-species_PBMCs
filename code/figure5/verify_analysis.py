"""Independent validation of GM-CSF filtering, paired effects and marker ranks."""
from pathlib import Path
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.stats import rankdata


def verify(work, checkpoint, source):
    work, checkpoint, source = map(Path, [work, checkpoint, source])
    table = pd.read_csv(work / "summary/population_changes.csv")
    donor = pd.read_csv(work / "summary/donor_population_scores_pSTAT5.csv")
    pairs = pd.read_csv(work / "summary/paired_donor_changes.csv")
    markers = pd.read_csv(work / "summary/T_cell_QBC_marker_trends.csv")
    for sp in ["human", "rhesus"]:
        old = pd.read_csv(source / f"{sp}_metadata.csv.gz")
        m = pd.read_csv(work / "inputs" / f"{sp}_metadata.csv.gz")
        expected = old.loc[old.condition.isin(["Basal", "GM-CSF"])]
        assert list(m.cell_id) == list(expected.cell_id)
        x = pd.read_csv(work / "inputs" / f"{sp}_matrix.csv.gz",
                        float_precision="round_trip").set_index("cell_id")
        oldx = pd.read_csv(source / f"{sp}_matrix.csv.gz",
                           float_precision="round_trip").set_index("cell_id")
        np.testing.assert_allclose(x, oldx.loc[m.cell_id], atol=1e-12, rtol=1e-12)
        scales = pd.read_csv(work / "inputs" / f"{sp}_signature_scale.csv")
        features = scales.feature.tolist()
        assert len(features) == 9 and not m.condition.eq("CD40L").any()
        recomputed = ((x[features].to_numpy() - scales["mean"].to_numpy()) /
                      scales.sd.to_numpy()).mean(axis=1)
        np.testing.assert_allclose(m.nine_marker_score, recomputed, atol=1e-12, rtol=1e-12)
        # Build a reference-cell table independently and compare donor summaries.
        check = m[["reference_id", "cytof_donor", "population", "condition"]].copy()
        check["nine_marker_score"] = recomputed
        check["pSTAT5"] = x["pSTAT5.Y694.-CyTOF"].to_numpy()
        ref = check.groupby(["reference_id", "cytof_donor", "population", "condition"]).mean(
            numeric_only=True).reset_index()
        for row in donor[donor.species == sp].itertuples():
            block = ref[(ref.cytof_donor == row.cytof_donor) &
                        (ref.population == row.population) &
                        (ref.condition == row.reference_group)]
            assert len(block) == row.n_unique_reference_cells
            for value in ["nine_marker_score", "pSTAT5"]:
                np.testing.assert_allclose(block[value].mean(), getattr(row, value),
                                           atol=1e-12, rtol=1e-12)
        for cond in ["Basal", "GM-CSF"]:
            cells = pd.read_csv(work / "graphs/tables" / f"{sp}_{cond}_cells.csv.gz")
            assert cells.condition.eq(cond).all()
            assert set(cells.cell_id) == set(m.loc[m.condition == cond, "cell_id"])
            np.testing.assert_array_equal(cells.root_cell,
                cells.nine_marker_score >= cells.nine_marker_score.quantile(.95))
        t = cells[cells.is_T & cells.finite_pseudotime]
        rank_pt = rankdata(t.pseudotime)
        saved = markers[markers.species == sp].sort_values("rank_within_species")
        assert len(saved) == 134 and saved.n_T_cells.eq(len(t)).all()
        for r in saved.itertuples():
            v = x.loc[t.cell_id, r.feature].to_numpy()
            rho = np.corrcoef(rankdata(v), rank_pt)[0, 1]
            np.testing.assert_allclose(rho, r.spearman_rho, atol=1e-12, rtol=1e-12)
        assert np.all(np.diff(saved.abs_rho) <= 1e-12)
    for r in table.itertuples():
        d = donor[(donor.species == r.species) & (donor.population == r.population)]
        c = d.pivot(index="cytof_donor", columns="reference_group",
                    values="n_unique_reference_cells").reindex(columns=["Basal", "GM-CSF"]).fillna(0)
        v = d.pivot(index="cytof_donor", columns="reference_group",
                    values=r.measure).reindex(columns=["Basal", "GM-CSF"])
        use = v.notna().all(axis=1) & c.ge(3).all(axis=1)
        z = v.loc[use, "GM-CSF"] - v.loc[use, "Basal"]
        assert len(z) == r.n_paired_donors
        assert r.eligible_for_ranking == (len(z) >= 5)
        for field, val in [("mean_delta", z.mean()), ("median_delta", z.median()),
                           ("donor_sd", z.std(ddof=1))]:
            np.testing.assert_allclose(val, getattr(r, field), atol=1e-12,
                                       rtol=1e-12, equal_nan=True)
        saved = pairs[(pairs.species == r.species) & (pairs.population == r.population) &
                      (pairs.measure == r.measure)].set_index("cytof_donor")
        assert set(saved.index) == set(z.index)
        np.testing.assert_allclose(saved.loc[z.index, "delta"], z, atol=1e-12, rtol=1e-12)
    for _, d in table[table.eligible_for_ranking].groupby(["species", "measure"]):
        np.testing.assert_array_equal(d.rank_absolute_change,
            d.mean_delta.abs().rank(method="min", ascending=False))
        np.testing.assert_array_equal(d.rank_signed_change,
            d.mean_delta.rank(method="min", ascending=False))
    copied = 0
    for directory in ["inputs", "graphs", "summary"]:
        sources = [p for p in (work / directory).rglob("*") if p.is_file()]
        targets = [p for p in (checkpoint / directory).rglob("*") if p.is_file()]
        assert {p.relative_to(work) for p in sources} == {
            p.relative_to(checkpoint) for p in targets}
        for p in sources:
            q = checkpoint / p.relative_to(work)
            assert p.stat().st_size == q.stat().st_size > 0
            assert hashlib.sha256(p.read_bytes()).digest() == hashlib.sha256(q.read_bytes()).digest()
            copied += 1
    result = {"status": "PASS", "CD40L_rows": 0,
              "score_scales_and_source_values": "PASS", "reference_donor_summaries": "PASS",
              "population_rankings": "PASS", "roots_and_all_134_marker_ranks": "PASS",
              "checkpoint_files_size_and_sha256_verified": copied}
    (work / "summary/analysis_validation.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
    return result
