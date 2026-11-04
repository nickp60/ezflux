#!/usr/bin/env python3
"""Match species abundances to the models of a (species-level) micom database.

Reads <sample>.annotated_abundances.csv (placeholder, sample_id, species, abundance) and writes
<prefix>.taxonomy.csv (micom Batch taxonomy) and <prefix>.annotation.tsv (models found per species).
"""
import argparse
import re
import sys

import pandas as pd
from micom.qiime_formats import load_qiime_manifest

RANKS = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]


def normalise(name):
    """Compare species names ignoring case, spaces/underscores and GTDB letter suffixes (Clostridium_A)."""
    name = re.sub(r"_[A-Z]{1,2}(?=\b|$)", "", str(name))
    return re.sub(r"[^a-z0-9]+", " ", name.lower()).strip()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--abundances", required=True)
    p.add_argument("--model-db", required=True)
    p.add_argument("--sample", required=True)
    p.add_argument("--prefix", required=True)
    args = p.parse_args()

    models = load_qiime_manifest(args.model_db).reset_index(drop=True)
    if models["summary_rank"].iloc[0] != "species":
        sys.exit(f"the model database is {models['summary_rank'].iloc[0]}-level; only species-level databases are supported")
    models["key"] = models["species"].map(normalise)

    abund = pd.read_csv(args.abundances, dtype={"sample_id": str})
    abund = abund[abund["abundance"] > 0].assign(key=lambda d: d["species"].map(normalise))
    abund["n_models"] = abund["key"].map(models["key"].value_counts()).fillna(0).astype(int)
    abund[["species", "abundance", "n_models"]].to_csv(f"{args.prefix}.annotation.tsv", sep="\t", index=False)

    found = abund.merge(models[RANKS + ["id", "key"]], on="key", suffixes=("_annotated", ""))
    if found.empty:
        sys.exit(f"none of the {len(abund)} species of sample {args.sample} are in the model database")
    # a species with several models shares its abundance equally; micom ids only allow [A-Za-z0-9_]
    found["abundance"] = found["abundance"] / found["n_models"]
    found["abundance"] = found["abundance"] / found["abundance"].sum()
    found["id"] = found["id"].str.replace(r"[^A-Za-z0-9_]+", "_", regex=True)
    found["sample_id"] = args.sample
    found[RANKS + ["id", "abundance", "sample_id"]].to_csv(f"{args.prefix}.taxonomy.csv", index=False)
    print(f"{args.sample}: {len(found)} models, {100 * abund.loc[abund.n_models > 0, 'abundance'].sum() / abund['abundance'].sum():.1f}% of abundance matched")


if __name__ == "__main__":
    main()
