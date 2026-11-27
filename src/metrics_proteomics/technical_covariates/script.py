"""Residual dependence on cell area. Thin wrapper over `sprotmetrics`."""

import argparse

import anndata as ad
import numpy as np
import sprotmetrics as spm

QC_COLUMN = "passed_QC"


def run(par: dict):
    layer = par.get("layer", "normalized")
    covariate = par.get("covariate", "area")
    n_pcs = int(par.get("n_pcs", 10))
    max_cells = int(par.get("max_cells", 5000))
    random_state = int(par.get("random_state", 0))

    adata = ad.read_h5ad(par["input"])
    scores = spm.covariate_dependence(
        adata,
        covariate,
        layer=layer,
        counts_layer="counts",
        n_pcs=n_pcs,
        max_cells=max_cells,
        random_state=random_state,
        qc_variants=True,
        qc_column=QC_COLUMN,
    )

    out = ad.AnnData(
        uns={
            "metric_ids": np.array(list(scores), dtype=object),
            "metric_values": np.array(list(scores.values()), dtype=float),
            "metric_group": "technical_covariates",
            "evaluated_layer": layer,
            "covariate": covariate,
            "covariate_present": covariate in adata.obs.columns,
            "counts_layer_present": "counts" in adata.layers,
            "n_pcs": n_pcs,
            "qc_column": QC_COLUMN,
        }
    )
    out.write_h5ad(par["output"], compression="gzip")


## VIASH START
par = {
    "input": "cell_by_protein_norm.h5ad",
    "layer": "normalized",
    "covariate": "area",
    "n_pcs": 10,
    "max_cells": 5000,
    "random_state": 0,
    "output": "metrics_technical_covariates.h5ad",
}
meta = {"name": "technical_covariates"}
## VIASH END


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--layer", default="normalized")
    ap.add_argument("--covariate", default="area")
    ap.add_argument("--n_pcs", type=int, default=10)
    ap.add_argument("--max_cells", type=int, default=5000)
    ap.add_argument("--random_state", type=int, default=0)
    ap.add_argument("--output", required=True)
    run(vars(ap.parse_args()))


if __name__ == "__main__":
    if "par" in globals():
        run(par)
    else:
        main()
