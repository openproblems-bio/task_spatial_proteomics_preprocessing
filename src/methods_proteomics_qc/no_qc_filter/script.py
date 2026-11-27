import anndata as ad
import numpy as np
import pandas as pd


def _get_counts_matrix(adata: ad.AnnData):
    return adata.layers["counts"] if "counts" in adata.layers else adata.X


def _get_area(adata: ad.AnnData):
    for col in ["area", "cell_area", "volume"]:
        if col in adata.obs.columns:
            return pd.to_numeric(adata.obs[col], errors="coerce").astype(float), col
    return pd.Series(np.nan, index=adata.obs_names, dtype=float), None


def _total_intensity(adata: ad.AnnData):
    X = _get_counts_matrix(adata)
    return pd.Series(np.asarray(X.sum(axis=1)).ravel().astype(float), index=adata.obs_names)


def run(par: dict):
    adata = ad.read_h5ad(par["input"])
    out_path = par["output"]

    area, area_col = _get_area(adata)
    total = _total_intensity(adata)
    passed = pd.Series(True, index=adata.obs_names)

    adata.obs["passed_QC"] = passed.values.astype(bool)
    adata.obs["total_intensity"] = total.values
    if area_col is None:
        adata.obs["area"] = area.values
    adata.obs["qc_area_ok"] = True
    adata.obs["qc_total_ok"] = True
    adata.obs["qc_filter_reason"] = "pass_no_qc"

    adata.uns["qc"] = {
        "method": "no_qc_filter",
        "description": "All cells are marked as passing QC.",
        "area_column_used": area_col if area_col is not None else "none",
        "n_cells": int(adata.n_obs),
        "n_passed": int(passed.sum()),
        "pass_rate": float(passed.mean()) if adata.n_obs > 0 else float("nan"),
        "median_total_intensity": float(np.nanmedian(total.values)) if adata.n_obs > 0 else float("nan"),
        "median_area": float(np.nanmedian(area.values)) if adata.n_obs > 0 else float("nan"),
    }

    adata.write_h5ad(out_path, compression="gzip")


if __name__ == "__main__":
    ## VIASH START
    par = {
        "input": None,
        "min_area": 50,
        "min_total": 0.0,
        "area_quantile": 0.05,
        "total_quantile": 0.05,
        "output": None,
    }
    meta = {}
    ## VIASH END
    run(par)
