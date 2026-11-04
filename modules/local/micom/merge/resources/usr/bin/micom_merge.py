#!/usr/bin/env python3
"""Merge the growth results of all samples and media, export tables and draw the plots.

Reads every *.zip (GrowthResults, one per sample and medium, sample_id = <sample>__<medium>) and
every *.build_manifest.csv in the working directory, so the files must be staged side by side.
Each plot is written as micom's interactive HTML and as a static *_mqc.png (picked up by MultiQC).
"""
import argparse
import glob
import sys
from functools import reduce
from itertools import cycle

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from micom.batch.results import load_results, save_results
from micom.measures import consumption_rates, production_rates
from micom.viz import plot_exchanges_per_taxon, plot_growth

MAX_TAXA_SHOWN = 60


def save_png(fig, path):
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def failed_plot(png, html, title, err):
    """Keep the pipeline going when a plot cannot be drawn (e.g. too few samples for UMAP)."""
    print(f"WARNING: {title} could not be drawn ({type(err).__name__}: {err})", file=sys.stderr)
    if html:
        with open(html, "w") as fh:
            fh.write(f"<html><body><p>{title} could not be drawn: {err}</p></body></html>")
    fig, ax = plt.subplots(figsize=(6, 1.5))
    ax.axis("off")
    ax.text(0.5, 0.5, f"{title} could not be drawn:\n{err}", ha="center", va="center", wrap=True)
    save_png(fig, png)


def niche_plot(growth, direction, html, png, n_neighbors):
    viz = plot_exchanges_per_taxon(growth, direction=direction, filename=html, n_neighbors=n_neighbors)
    umap = viz.data["reduced"]
    colors = dict(zip(sorted(umap.taxon.unique()), cycle(plt.get_cmap("tab20").colors)))
    fig, ax = plt.subplots(figsize=(8, 6))
    for taxon, rows in umap.groupby("taxon"):
        ax.scatter(rows["UMAP 1"], rows["UMAP 2"], s=28, alpha=0.8, color=colors[taxon], label=taxon)
    ax.set(xlabel="UMAP 1", ylabel="UMAP 2", title=f"Exchange niches: {direction} fluxes per taxon")
    if umap.taxon.nunique() <= 25:
        ax.legend(fontsize=7, loc="center left", bbox_to_anchor=(1.01, 0.5), frameon=False)
    save_png(fig, png)


def growth_plot(growth, html, png):
    rates = plot_growth(growth, filename=html).data["growth_rates"]
    order = rates.groupby("taxon").growth_rate.median().sort_values(ascending=False)
    shown = order.index[:MAX_TAXA_SHOWN]
    fig, ax = plt.subplots(figsize=(8, 0.3 * len(shown) + 2))
    rng = np.random.default_rng(0)
    for i, taxon in enumerate(shown):
        values = rates.loc[rates.taxon == taxon, "growth_rate"]
        ax.scatter(values, i + rng.uniform(-0.25, 0.25, len(values)), s=14, alpha=0.7)
    ax.set_yticks(range(len(shown)), shown, fontsize=7)
    ax.invert_yaxis()
    suffix = "" if len(order) <= len(shown) else f" (top {len(shown)} of {len(order)} taxa by median)"
    ax.set(xlabel="growth rate [1/h]", title="Taxa growth rates" + suffix)
    save_png(fig, png)


def found_fraction_plot(found, png):
    found = found.sort_values("found_abundance_fraction")
    pos = np.arange(len(found))
    fig, ax = plt.subplots(figsize=(8, 0.28 * len(found) + 2))
    ax.barh(pos + 0.2, found["found_abundance_fraction"], height=0.4, label="abundance found in model database")
    ax.barh(pos - 0.2, found["found_fraction"], height=0.4, label="taxa found in model database")
    ax.set_yticks(pos, found["sample_id"], fontsize=7)
    ax.axvline(0.5, color="grey", ls=":", lw=1)
    ax.set(xlim=(0, 1), xlabel="fraction", title="Share of each community matched to the model database")
    ax.legend(fontsize=7, loc="lower right")
    save_png(fig, png)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--prefix", required=True)
    p.add_argument("--n-neighbors", type=int, help="UMAP n_neighbors (default: max(2, min(8, n_samples // 2)))")
    args = p.parse_args()

    zips = sorted(glob.glob("*.zip"))
    growth = reduce(lambda a, b: a + b, (load_results(z) for z in zips))
    save_results(growth, f"{args.prefix}.growth.zip")
    growth.exchanges.to_csv(f"{args.prefix}.exchanges.csv", index=False)
    growth.growth_rates.to_csv(f"{args.prefix}.growth_rates.csv", index=False)
    production_rates(growth).to_csv(f"{args.prefix}.production_rates.csv")
    consumption_rates(growth).to_csv(f"{args.prefix}.consumption_rates.csv")

    # micom's default n_neighbors = min(8, n_samples // 2) is rejected by UMAP for < 4 samples
    n_samples = growth.exchanges.sample_id.nunique()
    n_neighbors = args.n_neighbors or max(2, min(8, n_samples // 2))
    for direction in ("import", "export"):
        html, png = f"{args.prefix}.{direction}_niche.html", f"{direction}_niche_mqc.png"
        try:
            niche_plot(growth, direction, html, png, n_neighbors)
        except Exception as err:
            failed_plot(png, html, f"The {direction} niche plot", err)
    try:
        growth_plot(growth, f"{args.prefix}.growth_rates.html", "growth_rates_mqc.png")
    except Exception as err:
        failed_plot("growth_rates_mqc.png", f"{args.prefix}.growth_rates.html", "The growth rate plot", err)

    # the same community is built once per medium: keep one row per sample
    found = pd.concat([pd.read_csv(f, dtype={"sample_id": str}) for f in sorted(glob.glob("*.build_manifest.csv"))])
    found = found.drop_duplicates("sample_id")
    columns = ["sample_id", "found_taxa", "total_taxa", "found_fraction", "found_abundance_fraction"]
    found[columns].to_csv(f"{args.prefix}.found_fraction.csv", index=False)
    try:
        found_fraction_plot(found, "found_fraction_mqc.png")
    except Exception as err:
        failed_plot("found_fraction_mqc.png", None, "The found fraction plot", err)
    print(f"merged {len(zips)} result sets, {n_samples} sample ids")


if __name__ == "__main__":
    main()
