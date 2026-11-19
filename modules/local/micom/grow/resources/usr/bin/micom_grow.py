#!/usr/bin/env python3
"""Grow one community (one sample on one medium) with micom's Batch workflow.

Writes <out>.zip (GrowthResults, sample_id relabelled to --label), <out>.build_manifest.csv and the built
community models in --models-dir.
Exits with 1 when the sample could not be solved so Nextflow can retry.

With --collapse-rare-threshold the rare taxa are pooled into one meta-organism called "rare" before the community is built:
taxa below --cutoff are dropped, taxa between --cutoff and the threshold are collapsed into "rare" (the union of their
reactions, so the pathways they contribute stay available; biomass reaction = their abundance-weighted average), and taxa
at or above the threshold stay full organisms. --collapse-rare-abundance sets the abundance given to "rare": the sum of the
abundances of the collapsed taxa, their median, or the threshold itself (default).
"""
import argparse
import os
import re
import sys
import shutil
from pathlib import Path
from zipfile import ZipFile

import micom.community as mc
import micom.util as mu
import pandas as pd
from cobra import Metabolite, Reaction
from cobra.io import save_json_model
from cobra.util.solver import interface_to_str, linear_reaction_coefficients
from micom import Batch
from micom.batch import Configuration
from micom.batch.results import save_results
from micom.constants import RANKS
from micom.db import load_manifest, load_zip_model_db
from micom.qiime_formats import load_qiime_manifest, load_qiime_medium, metadata
from micom.taxonomy import unify_rank_prefixes
from micom.util import load_model

RARE_ID = "rare"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--taxonomy", required=True, help="micom taxonomy CSV")
    p.add_argument("--medium", required=True, help="medium: q2-micom artifact (.qza) or CSV/TSV with columns reaction, flux")
    p.add_argument("--model-db", required=True, help="micom model database")
    p.add_argument("--label", required=True, help="sample_id used in the results, e.g. <sample>__<medium>")
    p.add_argument("--out", required=True, help="output prefix")
    p.add_argument("--models-dir", default="models", help="folder for the built community models")
    p.add_argument("--solver", default="gurobi")
    p.add_argument("--threads", type=int, default=1)
    p.add_argument("--tradeoff", type=float, default=0.5)
    p.add_argument("--cutoff", type=float, default=0.001)
    p.add_argument("--tolerance", type=float, default=1e-6)
    p.add_argument("--no-coupling", action="store_true", help="switch off resource coupling constraints")
    p.add_argument("--presolve", choices=["true", "false", "auto"], help="override the solver presolve")
    p.add_argument("--numeric-focus", type=int, choices=[0, 1, 2, 3], help="Gurobi NumericFocus")
    p.add_argument("--bar-conv-tol", type=float, help="Gurobi BarConvTol (micom sets 1e-9)")
    p.add_argument("--bar-iter-limit", type=int, help="Gurobi BarIterLimit (micom sets 2000)")
    p.add_argument("--solver-log", action="store_true", help="let Gurobi write its own solver log (very verbose) to the task's stdout")
    p.add_argument("--collapse-rare-threshold", type=float, help="collapse the taxa with a relative abundance between --cutoff and this value into one meta-organism 'rare'; must be larger than --cutoff")
    p.add_argument("--collapse-rare-abundance", choices=["sum", "median", "threshold"], default="threshold",
                   help="abundance of the collapsed meta-organism: sum or median of the abundances of the collapsed taxa, or the --collapse-rare-threshold value itself")
    p.add_argument("--flux_method", type=str, choices=['minimal imports', 'pFBA', 'none'] , help="MICOM flux method", default='pFBA')
    return p.parse_args()


def load_medium(path):
    """A q2-micom medium (.qza), or a CSV/TSV with the same table: columns `reaction` and `flux` (mmol/gDW/h)."""
    if Path(path).suffix == ".qza":
        return load_qiime_medium(path)
    medium = pd.read_csv(path, sep=None, engine="python")   # comma or tab
    missing = {"reaction", "flux"} - set(medium.columns)
    if missing:
        sys.exit(f"{path}: medium table needs the columns reaction and flux, missing {sorted(missing)}")
    if medium.reaction.isna().any() or medium.reaction.duplicated().any():
        sys.exit(f"{path}: reaction ids must be unique and not empty")
    medium["flux"] = pd.to_numeric(medium.flux, errors="coerce")
    if medium.flux.isna().any() or (medium.flux <= 0).any():
        sys.exit(f"{path}: flux must be a positive number for every reaction")
    medium.index = medium.reaction
    return medium


# micom hard-codes its Gurobi settings in adjust_solver_config() and Batch exposes none of them.
# Batch workers are *spawned* processes that re-import this file (with the same sys.argv), so the
# override must be applied at import time, not inside main().
ARGS = parse_args()
PRESOLVE = {"true": True, "false": False, "auto": "auto", None: None}[ARGS.presolve]
_micom_adjust = mu.adjust_solver_config


GUROBI_PARAMS = ["Method", "BarConvTol", "BarIterLimit", "NumericFocus", "Presolve", "ScaleFlag", "CrossoverBasis", "Threads"]
_LOGGED = set()


def _log_settings(solver):
    """Print the Gurobi parameters that are set on the model, once per process and distinct set of values."""
    params = solver.problem.Params
    values = tuple((name, getattr(params, name, None)) for name in GUROBI_PARAMS)
    if values in _LOGGED:
        return
    _LOGGED.add(values)
    shown = ", ".join(f"{name}={value}" for name, value in values)
    print(f"[micom_grow] Gurobi settings in effect (pid {os.getpid()}): {shown}", flush=True)


def _adjust(solver):
    _micom_adjust(solver)
    if PRESOLVE is not None:
        solver.configuration.presolve = PRESOLVE
    if interface_to_str(solver.interface) == "gurobi":
        if ARGS.numeric_focus is not None:
            solver.problem.Params.NumericFocus = ARGS.numeric_focus
        if ARGS.bar_conv_tol is not None:
            solver.problem.Params.BarConvTol = ARGS.bar_conv_tol
        if ARGS.bar_iter_limit is not None:
            solver.problem.Params.BarIterLimit = ARGS.bar_iter_limit
        if ARGS.solver_log:
            solver.problem.Params.LogToConsole = 1   # micom switches the console log off
        _log_settings(solver)


mu.adjust_solver_config = _adjust
mc.adjust_solver_config = _adjust


def db_manifest(model_db, folder):
    """The model database manifest with a `file` column of paths, and a function that makes sure the files exist in `folder`."""
    path = Path(model_db)
    if path.suffix == ".qza":   # only the models that are needed are extracted
        uuid = metadata(path)["uuid"]
        manifest = load_qiime_manifest(path)
        manifest["file"] = [os.path.join(folder, uuid, "data", f) for f in manifest.file]

        def fetch(files):
            with ZipFile(path) as zf:
                for f in files:
                    zf.extract(os.path.relpath(f, folder), folder)

        return manifest, fetch
    manifest = load_zip_model_db(path, folder) if path.suffix == ".zip" else load_manifest(path)
    return manifest, lambda files: None


def join_weighted(files, weights, model_id):
    """Union of the reactions of several models with one biomass reaction: the abundance-weighted average of their biomass reactions.

    The models must share one ID system (as the models of a database do). Where several models have the same reaction,
    the first one (the most abundant) is kept. The members' own biomass reactions are replaced by the combined one.
    """
    models = [load_model(f) for f in files]
    model = models[0]
    stoichiometry, templates = {}, {}
    for weight, member in zip(weights, models):
        for rxn, coef in linear_reaction_coefficients(member, member.reactions).items():
            for met, c in rxn.metabolites.items():
                stoichiometry[met.id] = stoichiometry.get(met.id, 0.0) + weight * coef * c
                templates.setdefault(met.id, met)
    own_biomass = [r.id for r in linear_reaction_coefficients(model, model.reactions)]
    have = {r.id for r in model.reactions}
    for member in models[1:]:
        biomass_ids = {r.id for r in linear_reaction_coefficients(member, member.reactions)}
        new = [r for r in member.reactions if r.id not in have and r.id not in biomass_ids]
        model.add_reactions(new)
        have.update(r.id for r in new)
    model.remove_reactions([model.reactions.get_by_id(i) for i in own_biomass], remove_orphans=False)
    biomass = Reaction(id=f"{model_id}_biomass", name="abundance-weighted biomass of the collapsed taxa",
                       subsystem="biomass production", lower_bound=0, upper_bound=1000)
    model.add_reactions([biomass])
    biomass.add_metabolites({
        (model.metabolites.get_by_id(i) if i in model.metabolites
         else Metabolite(i, formula=templates[i].formula, name=templates[i].name, compartment=templates[i].compartment,
                         charge=templates[i].charge)): c
        for i, c in stoichiometry.items() if c != 0
    })
    model.objective = biomass
    model.id = model_id
    return model


def collapse_rare(taxonomy, model_db, cutoff, threshold, folder, abundance_mode="threshold"):
    """Pool the rare taxa of the (single) sample into one meta-organism.

    Follows micom's own order: abundances are made relative, taxa <= cutoff are dropped, the rest is matched to the
    model database. Of the matched taxa those below the threshold are joined into `rare`, whose abundance is the sum or
    the median of the abundances of those taxa, or the threshold itself (abundance_mode). As always the abundances of the
    community are made relative again when the community is built.
    Returns (taxonomy with a `file` column, build metrics, table of the collapsed taxa), or None when fewer than two
    taxa are below the threshold (nothing to collapse; the normal model database route is used then).
    """
    tax = taxonomy.drop(columns=["file"], errors="ignore").copy()
    if "abundance" not in tax.columns:
        tax["abundance"] = 1
    tax["abundance"] = tax.abundance / tax.abundance.sum()
    tax = tax[tax.abundance > cutoff]
    manifest, fetch = db_manifest(model_db, folder)
    rank = manifest["summary_rank"][0]
    if rank not in tax.columns:
        sys.exit(f"Missing the column `{rank}` from the taxonomy.")
    if "id" not in tax.columns:
        tax["id"] = tax[rank]
    keep = [r for r in RANKS[0:RANKS.index(rank) + 1] if r in tax.columns and r in manifest.columns]
    manifest = manifest[keep + ["file"]]
    tax = unify_rank_prefixes(tax, manifest)
    merged = pd.merge(tax, manifest, on=keep)
    metrics = {
        "found_taxa": merged.shape[0],
        "total_taxa": tax.shape[0],
        "found_fraction": merged.shape[0] / tax.shape[0] if tax.shape[0] else 0.0,
        "found_abundance_fraction": merged.abundance.sum(),
    }
    merged["id"] = merged.id.str.replace(r"[^A-Za-z0-9_]+", "_", regex=True)
    merged = merged.sort_values("abundance", ascending=False)
    rare = merged[merged.abundance < threshold]
    full = merged[merged.abundance >= threshold]
    print(f"Taxa matched in the model database: {len(merged)}; full organisms (abundance >= {threshold}): {len(full)}; "
          f"rare (cutoff < abundance < {threshold}): {len(rare)}", flush=True)
    if len(rare) < 2:
        print("Fewer than two rare taxa: nothing to collapse.", flush=True)
        return None
    rare_id = RARE_ID
    while rare_id in set(full.id):
        rare_id += "_"
    weights = (rare.abundance / rare.abundance.sum()).tolist()
    rare_abundance = {"sum": rare.abundance.sum(), "median": rare.abundance.median(), "threshold": threshold}[abundance_mode]
    rare_fraction = rare_abundance / (full.abundance.sum() + rare_abundance)   # share of the finished community
    fetch(list(full.file) + list(rare.file))
    rare_file = os.path.join(folder, f"{rare_id}.json")
    model = join_weighted(list(rare.file), weights, rare_id)
    save_json_model(model, rare_file)
    print(f"Collapsed {len(rare)} taxa ({100 * rare.abundance.sum() / merged.abundance.sum():.1f}% of the matched abundance) into "
          f"'{rare_id}': {len(model.reactions)} reactions, {len(model.metabolites)} metabolites", flush=True)
    print(f"Abundance of '{rare_id}' ({abundance_mode} of the collapsed taxa): {rare_abundance:.5g}, "
          f"{100 * rare_fraction:.1f}% of the community", flush=True)
    rare_row = pd.DataFrame({"id": [rare_id], "abundance": [rare_abundance], "file": [rare_file]})
    rows = pd.concat([full[["id", "abundance", "file"]], rare_row], ignore_index=True)
    rows.insert(0, "sample_id", taxonomy.sample_id.iloc[0])
    collapsed = pd.DataFrame({"id": rare.id.values, "abundance": rare.abundance.values, "weight_in_rare": weights})
    metrics.update({"n_full_taxa": len(full), "n_collapsed_taxa": len(rare),
                    "collapsed_abundance_fraction": rare.abundance.sum() / merged.abundance.sum(),
                    "rare_abundance_mode": abundance_mode, "rare_abundance_fraction": rare_fraction})
    return rows, metrics, collapsed


def main():
    config = Configuration()
    config.solver = ARGS.solver
    config.threads = ARGS.threads
    config.tolerance = ARGS.tolerance
    config.progress = False
    config.simulation.tradeoff = ARGS.tradeoff
    config.simulation.flux_method = ARGS.flux_method
    config.coupling.enabled = not ARGS.no_coupling
    config.build.cutoff = ARGS.cutoff
    config.build.force_rebuild = True
    print(config)

    taxonomy = pd.read_csv(ARGS.taxonomy, dtype={"sample_id": str})
    if taxonomy.sample_id.nunique() != 1:
        sys.exit(f"{ARGS.taxonomy} must hold exactly one sample, found {taxonomy.sample_id.nunique()}")

    model_db, metrics, collapsed = ARGS.model_db, None, None
    workdir = Path(f"{ARGS.out}.collapse_work").resolve()   # in the task's work directory, not in the system temp
    workdir.mkdir(exist_ok=True)
    try:
        if ARGS.collapse_rare_threshold is not None:
            if not ARGS.collapse_rare_threshold > ARGS.cutoff:
                sys.exit(f"--collapse-rare-threshold ({ARGS.collapse_rare_threshold}) must be larger than --cutoff ({ARGS.cutoff})")
            result = collapse_rare(taxonomy, ARGS.model_db, ARGS.cutoff, ARGS.collapse_rare_threshold, str(workdir), ARGS.collapse_rare_abundance)
            if result is not None:
                taxonomy, metrics, collapsed = result
                model_db = None            # the taxonomy carries the model files now
                config.dbs.microbial = None
        run(config, taxonomy, model_db, metrics, collapsed)
    finally:
        shutil.rmtree(workdir, ignore_errors=True)   # the extracted models; the built community pickle holds what is needed


def run(config, taxonomy, model_db, metrics, collapsed):
    batch = Batch(
        taxonomy=taxonomy,
        medium=load_medium(ARGS.medium),
        model_db=model_db,
        config=config,
        name=ARGS.out,
    )
    batch.build(ARGS.models_dir)
    if metrics is not None:   # same columns as the model database route, plus what was collapsed
        for key, value in metrics.items():
            batch.build_manifest[key] = value
        batch.build_manifest.to_csv(Path(ARGS.models_dir) / "manifest.csv", index=False)
        collapsed.to_csv(f"{ARGS.out}.collapsed_taxa.csv", index=False)
    batch.build_manifest.to_csv(f"{ARGS.out}.build_manifest.csv", index=False)
    print("Build manifest:")
    print(batch.build_manifest.to_string(index=False), flush=True)   # before the (long) growth simulation
    results = batch.grow()   # raises if the sample cannot be solved

    results.growth_rates["sample_id"] = ARGS.label
    results.exchanges["sample_id"] = ARGS.label
    save_results(results, f"{ARGS.out}.zip")


if __name__ == "__main__":
    main()
