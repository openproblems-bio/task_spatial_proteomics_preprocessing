"""Raw-versus-normalized diagnostics. Thin wrapper over `sprotmetrics`."""

import argparse

import anndata as ad
import numpy as np
import sprotmetrics as spm

QC_COLUMN = "passed_QC"


def run(par: dict):
    raw_layer = par.get("raw_layer", "counts")
    normalized_layer = par.get("normalized_layer", "normalized")
    max_cells = int(par.get("max_cells", 5000))
    random_state = int(par.get("random_state", 0))

    adata = ad.read_h5ad(par["input"])
    scores = spm.normalization_diagnostics(
        adata,
        raw_layer=raw_layer,
        normalized_layer=normalized_layer,
        max_cells=max_cells,
        random_state=random_state,
        qc_variants=True,
        qc_column=QC_COLUMN,
    )

    out = ad.AnnData(
        uns={
            "metric_ids": np.array(list(scores), dtype=object),
            "metric_values": np.array(list(scores.values()), dtype=float),
            "metric_group": "normalization_diagnostics",
            "raw_layer": raw_layer,
            "normalized_layer": normalized_layer,
            "max_cells": max_cells,
            "random_state": random_state,
            "qc_column": QC_COLUMN,
        }
    )
    out.write_h5ad(par["output"], compression="gzip")


## VIASH START
par = {
    "input": "cell_by_protein_norm.h5ad",
    "raw_layer": "counts",
    "normalized_layer": "normalized",
    "max_cells": 5000,
    "random_state": 0,
    "output": "metrics_normalization_diagnostics.h5ad",
}
meta = {"name": "normalization_diagnostics"}
## VIASH END


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--raw_layer", default="counts")
    ap.add_argument("--normalized_layer", default="normalized")
    ap.add_argument("--max_cells", type=int, default=5000)
    ap.add_argument("--random_state", type=int, default=0)
    ap.add_argument("--output", required=True)
    run(vars(ap.parse_args()))


if __name__ == "__main__":
    if "par" in globals():
        run(par)
    else:
        main()
