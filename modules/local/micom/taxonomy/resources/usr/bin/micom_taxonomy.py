#!/usr/bin/env python3
"""Match the taxa of one sample to the models of a (species-level) micom database.

Two kinds of input, both ending in the same micom Batch taxonomy table:
  --source sylph    the raw `sylph profile` table. Every genome (Genome_file) is looked up in the model
                    metadata by its path, which gives the species name of that genome; its Taxonomic_abundance
                    goes to the model of that species. Genomes are never merged by name beforehand.
  --source species  a table with the columns species and abundance (the ASV/BLAST route).

A species name is matched to the `species` column of the model database exactly, so "[Eubacterium] rectale"
and "Eubacterium rectale" stay two different models. Only when there is no exact match are brackets and
underscores ignored. Writes <prefix>.taxonomy.csv and <prefix>.annotation.tsv (what happened to every input row).
"""
import argparse
import re
import sys

import pandas as pd
from micom.qiime_formats import load_qiime_manifest

RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]


def loose(name):
    """Fallback key: no brackets, underscores read as spaces."""
    return re.sub(r"\s+", " ", re.sub(r"[\[\]]", "", str(name)).replace("_", " ")).strip()


def read_metadata(path):
    """Genome_file -> species of the bacterial genomes in the model metadata (no header, tab separated)."""
    meta = pd.read_csv(path, sep="\t", header=None, dtype=str)
    if meta.shape[1] != 2:
        sys.exit(f"{path}: the model metadata needs two tab-separated columns, Genome_file and clade_name")
    meta.columns = ["Genome_file", "clade_name"]
    if meta.Genome_file.duplicated().any():
        sys.exit(f"{path}: Genome_file is listed more than once: {meta.Genome_file[meta.Genome_file.duplicated()].unique()[:5].tolist()}")
    meta["species"] = meta.clade_name.str.extract(r"s__(.*)$", expand=False).str.strip()
    meta["bacterial"] = meta.clade_name.str.startswith("d__Bacteria")
    meta["base"] = meta.Genome_file.map(lambda f: f.rsplit("/", 1)[-1])
    return meta


def sylph_to_species(profile_path, meta):
    """One row per sylph genome: the species of its Genome_file, looked up by path (then by unique file name)."""
    prof = pd.read_csv(profile_path, sep="\t", dtype={"Genome_file": str})
    if not {"Genome_file", "Taxonomic_abundance"} <= set(prof.columns):
        sys.exit(f"{profile_path} is not a sylph profile (needs the columns Genome_file and Taxonomic_abundance)")
    by_path = meta.set_index("Genome_file")
    unique_base = meta[~meta.base.duplicated(keep=False)].set_index("base")
    rows = []
    for genome, abundance in zip(prof.Genome_file, prof.Taxonomic_abundance):
        if genome in by_path.index:
            hit, how = by_path.loc[genome], "path"
        elif genome.rsplit("/", 1)[-1] in unique_base.index:
            hit, how = unique_base.loc[genome.rsplit("/", 1)[-1]], "file name"
        else:
            rows.append((genome, None, abundance, "not in model metadata"))
            continue
        if not hit.bacterial:
            rows.append((genome, hit.species, abundance, "not bacterial"))
        elif not isinstance(hit.species, str) or not hit.species:
            rows.append((genome, None, abundance, "no species in clade_name"))
        else:
            rows.append((genome, hit.species, abundance, f"metadata ({how})"))
    return pd.DataFrame(rows, columns=["entry", "species", "abundance", "status"])


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--abundances", required=True)
    p.add_argument("--source", choices=["sylph", "species"], required=True)
    p.add_argument("--model-metadata", help="Genome_file <tab> clade_name (needed with --source sylph)")
    p.add_argument("--model-db", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--prefix", required=True)
    args = p.parse_args()

    models = load_qiime_manifest(args.model_db).reset_index(drop=True)
    if models["summary_rank"].iloc[0] != "species":
        sys.exit(f"the model database is {models['summary_rank'].iloc[0]}-level; only species-level databases are supported")

    if args.source == "sylph":
        if not args.model_metadata:
            sys.exit("--source sylph needs --model-metadata")
        table = sylph_to_species(args.abundances, read_metadata(args.model_metadata))
    else:
        table = pd.read_csv(args.abundances, dtype={"sample_id": str})
        table = table.assign(entry=table["species"], status="species table")[["entry", "species", "abundance", "status"]]
    table = table[table["abundance"] > 0].reset_index(drop=True)

    # species name -> rows of the model database: exact first, then ignoring brackets and underscores
    exact = models.groupby("species").groups
    models["loose"] = models["species"].map(loose)
    fuzzy = models.groupby("loose").groups
    hits, match = [], []
    for species, status in zip(table["species"], table["status"]):
        if species is None or status == "not bacterial":
            hits.append([]); match.append(status); continue
        if species in exact:
            hits.append(list(exact[species])); match.append("exact")
        elif loose(species) in fuzzy:
            hits.append(list(fuzzy[loose(species)])); match.append("ignoring brackets")
        else:
            hits.append([]); match.append("no model")
    table["n_models"] = [len(h) for h in hits]
    table["match"] = match
    table[["entry", "species", "abundance", "n_models", "match"]].to_csv(f"{args.prefix}.annotation.tsv", sep="\t", index=False)

    # every input row shares its abundance equally over the models of its species; models are then summed
    found = pd.DataFrame(
        [(m, a / len(h)) for h, a in zip(hits, table["abundance"]) for m in h], columns=["model", "abundance"]
    ).groupby("model", as_index=False).sum()
    if found.empty:
        sys.exit(f"none of the {len(table)} taxa of sample {args.sample} are in the model database")
    found["abundance"] = found["abundance"] / found["abundance"].sum()
    out = models.loc[found["model"], RANKS + ["id"]].reset_index(drop=True).assign(abundance=found["abundance"].values)
    # micom ids only allow [A-Za-z0-9_]; one character each, so "[Eubacterium] rectale" and "Eubacterium rectale" stay apart
    out["id"] = out["id"].str.replace(r"[^A-Za-z0-9_]", "_", regex=True)
    dup = out["id"][out["id"].duplicated()].tolist()
    if dup:
        sys.exit(f"model ids are not unique after making them safe for micom: {dup}")
    out["sample_id"] = args.sample
    out.to_csv(f"{args.prefix}.taxonomy.csv", index=False)
    matched = table.loc[table.n_models > 0, "abundance"].sum() / table["abundance"].sum()
    print(f"{args.sample}: {len(out)} models, {100 * matched:.1f}% of abundance matched")


if __name__ == "__main__":
    main()
