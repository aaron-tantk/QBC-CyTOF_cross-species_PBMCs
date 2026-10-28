"""Attached-style A–H Figure 5 using only Basal and GM-CSF references."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, Normalize
from matplotlib.lines import Line2D
from adjustText import adjust_text

plt.rcParams.update({"font.family": ["Liberation Sans", "Arimo", "DejaVu Sans"],
                     "font.size": 15, "axes.titlesize": 18,
                     "axes.spines.top": False, "axes.spines.right": False})
NAMES = {"CD4_T_cells": "CD4 T", "CD8_T_cells": "CD8 T",
         "Activated_T_cells": "Activated T", "Classical_monocytes": "Monocytes",
         "Dendritic_cells": "DC", "Plasmacytoid_Dendritic_cells": "pDC",
         "NK_cells": "NK", "B_cells": "B", "Activated_dendritic_cells": "Activated DC",
         "NK-NKT_cells": "NK/NKT"}
TYPES = ["CD4_T_cells", "CD8_T_cells", "Activated_T_cells"]
TCOLORS = ["#0072B2", "#E69F00", "#009E73"]


def labels(ax, d):
    texts = []
    for name, g in d.groupby("population"):
        h, xe, ye = np.histogram2d(g.UMAP1, g.UMAP2, bins=12)
        i, j = np.unravel_index(np.argmax(h), h.shape)
        local = g[g.UMAP1.between(xe[i], xe[i+1]) & g.UMAP2.between(ye[j], ye[j+1])]
        xy = local[["UMAP1", "UMAP2"]].median()
        texts.append(ax.text(*xy, NAMES.get(name, name.replace("_", " ")),
            ha="center", va="center", fontsize=12.5,
            bbox={"boxstyle": "round,pad=.15", "fc": "white", "ec": ".6",
                  "alpha": .9, "lw": .5}))
    np.random.seed(42)
    adjust_text(texts, ax=ax, force_text=(.6, .8), expand=(1.05, 1.2),
                ensure_inside_axes=True, time_lim=1.5,
                arrowprops={"arrowstyle": "-", "color": ".3", "lw": .6})


def scatter(ax, d, field, title, cmap, norm=None, colorbar=True):
    finite = np.isfinite(d[field])
    ax.scatter(d.loc[~finite, "UMAP1"], d.loc[~finite, "UMAP2"],
               s=3, c="#C8C8C8", rasterized=True)
    q = d.loc[finite].sample(frac=1, random_state=42)
    im = ax.scatter(q.UMAP1, q.UMAP2, c=q[field], s=3, cmap=cmap,
                    norm=norm, rasterized=True)
    ax.set(title=title, xlabel="UMAP 1", ylabel="UMAP 2", xticks=[], yticks=[])
    ax.margins(.12)
    labels(ax, d)
    if colorbar:
        cb = ax.figure.colorbar(im, ax=ax, fraction=.035, pad=.025)
        cb.set_label("Pseudotime", fontsize=13)
        cb.ax.tick_params(labelsize=12)
    return im


def heatmap(fig, slot, work, species, panel):
    d = pd.read_csv(work / "graphs/tables" / f"{species}_GM-CSF_cells.csv.gz")
    d = d[d.is_T & d.finite_pseudotime].sort_values("pseudotime", kind="stable")
    ranks = pd.read_csv(work / "summary/T_cell_QBC_marker_trends.csv")
    top = ranks[ranks.species == species].head(20).feature.tolist()
    x = pd.read_csv(work / "inputs" / f"{species}_matrix.csv.gz",
                    float_precision="round_trip").set_index("cell_id")
    z = x.loc[d.cell_id, top].to_numpy().T
    z = (z - z.mean(axis=1, keepdims=True)) / z.std(axis=1, ddof=1, keepdims=True)
    grid = slot.subgridspec(3, 2, height_ratios=[.13, .13, 4.4],
                           width_ratios=[1, .04], hspace=.07, wspace=.12)
    pt = fig.add_subplot(grid[0, 0])
    subset = fig.add_subplot(grid[1, 0])
    ax = fig.add_subplot(grid[2, 0])
    pt.imshow(d.pseudotime.to_numpy()[None, :], aspect="auto", cmap="viridis")
    subset.imshow(d.population.map({k:i for i,k in enumerate(TYPES)}).to_numpy()[None, :],
                  aspect="auto", cmap=ListedColormap(TCOLORS), vmin=0, vmax=2)
    for track, title in [(pt, "Pseudotime"), (subset, "T subset")]:
        track.set(xticks=[], yticks=[])
        track.set_ylabel(title, rotation=0, ha="right", va="center", fontsize=12, labelpad=6)
        for spine in track.spines.values():
            spine.set_visible(False)
    im = ax.imshow(z, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5,
                   interpolation="nearest")
    ax.set(yticks=range(len(top)), yticklabels=[c.removesuffix("-QBC") for c in top],
           xticks=[], xlabel="T cells ordered by pseudotime")
    ax.tick_params(axis="y", labelsize=13)
    cb = fig.colorbar(im, cax=fig.add_subplot(grid[2, 1]))
    cb.set_label("Within-marker z", fontsize=12)
    cb.ax.tick_params(labelsize=12)
    pt.set_title(f"{panel}  {species.title()} T cells — GM-CSF", loc="left", pad=20)


def plot(work, output):
    work, output = Path(work), Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(14, 21))
    outer = fig.add_gridspec(5, 1, height_ratios=[1, 1, 1, 1, 1.2], hspace=.42,
                            top=.97, bottom=.05, left=.09, right=.95)
    cache = {(s, c): pd.read_csv(work / "graphs/tables" / f"{s}_{c}_cells.csv.gz")
             for s in ["human", "rhesus"] for c in ["Basal", "GM-CSF"]}
    for row, sp, panels in [(0, "human", ["A", "B"]), (1, "rhesus", ["C", "D"])]:
        grid = outer[row].subgridspec(1, 2, wspace=.3)
        for j, cond in enumerate(["Basal", "GM-CSF"]):
            scatter(fig.add_subplot(grid[0,j]), cache[(sp,cond)], "pseudotime",
                    f"{panels[j]}  {sp.title()} — {cond}", "viridis")
    display = {}
    for row, sp, panel in [(2, "human", "E"), (3, "rhesus", "F")]:
        both = pd.concat([cache[(sp,c)] for c in ["Basal", "GM-CSF"]])
        low, high = both.nine_marker_score.quantile([.01, .99])
        display[sp] = {"score_color_min": float(low), "score_color_max": float(high),
                       "limits": "species_pooled_Basal_GM-CSF_1st_99th_percentiles"}
        grid = outer[row].subgridspec(1, 3, width_ratios=[1, 1, .035], wspace=.35)
        ims = []
        for j, cond in enumerate(["Basal", "GM-CSF"]):
            title = f"{panel}  {sp.title()} — {cond}" if j == 0 else f"{sp.title()} — {cond}"
            ims.append(scatter(fig.add_subplot(grid[0,j]), cache[(sp,cond)], "nine_marker_score",
                               title, "magma", Normalize(low, high), colorbar=False))
        cb = fig.colorbar(ims[0], cax=fig.add_subplot(grid[0,2]))
        cb.set_label("Nine-marker score", fontsize=13)
    grid = outer[4].subgridspec(1, 2, wspace=.75)
    heatmap(fig, grid[0,0], work, "human", "G")
    heatmap(fig, grid[0,1], work, "rhesus", "H")
    fig.legend(handles=[Line2D([0],[0],color=c,lw=6) for c in TCOLORS],
               labels=["CD4 T","CD8 T","Activated T"], loc="lower center",
               ncol=3, bbox_to_anchor=(.5,.005), frameon=False)
    fig.savefig(output, dpi=300, bbox_inches="tight", pad_inches=.12)
    plt.close(fig)
    import json
    (work / "summary/figure_display_scales.json").write_text(json.dumps(display, indent=2))
    print(output)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--work", required=True)
    ap.add_argument("--output", required=True)
    a = ap.parse_args()
    plot(a.work, a.output)
