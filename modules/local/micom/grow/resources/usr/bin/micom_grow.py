#!/usr/bin/env python3
"""Grow one community (one sample on one medium) with micom's Batch workflow.

Writes <out>.zip (GrowthResults, sample_id relabelled to --label), <out>.build_manifest.csv and the built
community models in --models-dir.
Exits with 1 when the sample could not be solved so Nextflow can retry.
"""
import argparse
import os
import sys
from pathlib import Path

import micom.community as mc
import micom.util as mu
import pandas as pd
from cobra.util.solver import interface_to_str
from micom import Batch
from micom.batch import Configuration
from micom.batch.results import save_results
from micom.qiime_formats import load_qiime_medium


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

    batch = Batch(
        taxonomy=taxonomy,
        medium=load_medium(ARGS.medium),
        model_db=ARGS.model_db,
        config=config,
        name=ARGS.out,
    )
    batch.build(ARGS.models_dir)
    batch.build_manifest.to_csv(f"{ARGS.out}.build_manifest.csv", index=False)
    print("Build manifest:")
    print(batch.build_manifest.to_string(index=False), flush=True)   # before the (long) growth simulation
    results = batch.grow()   # raises if the sample cannot be solved

    results.growth_rates["sample_id"] = ARGS.label
    results.exchanges["sample_id"] = ARGS.label
    save_results(results, f"{ARGS.out}.zip")


if __name__ == "__main__":
    main()
