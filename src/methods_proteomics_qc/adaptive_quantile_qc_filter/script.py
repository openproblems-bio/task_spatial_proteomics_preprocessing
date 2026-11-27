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


def _safe_quantile(values, q):
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return float("nan")
    return float(np.quantile(arr, q))


def _failure_reason(area_ok: pd.Series, total_ok: pd.Series):
    reasons = []
    for a_ok, t_ok in zip(area_ok.values, total_ok.values):
        if bool(a_ok) and bool(t_ok):
            reasons.append("pass")
        elif (not bool(a_ok)) and (not bool(t_ok)):
            reasons.append("low_area_quantile_and_low_total_quantile")
        elif not bool(a_ok):
            reasons.append("low_area_quantile")
        else:
            reasons.append("low_total_quantile")
    return reasons


def run(par: dict):
    adata = ad.read_h5ad(par["input"])
    out_path = par["output"]
    area_quantile = float(par.get("area_quantile", 0.05))
    total_quantile = float(par.get("total_quantile", 0.05))

    if not (0 <= area_quantile < 1):
        raise ValueError("area_quantile must be in [0, 1).")
    if not (0 <= total_quantile < 1):
        raise ValueError("total_quantile must be in [0, 1).")

    area, area_col = _get_area(adata)
    total = _total_intensity(adata)

    if area_col is None:
        area_threshold = float("nan")
        area_ok = pd.Series(True, index=adata.obs_names)
    else:
        area_threshold = _safe_quantile(area.values, area_quantile)
        area_ok = (area >= area_threshold).fillna(False)

    total_threshold = _safe_quantile(total.values, total_quantile)
    total_ok = (total >= total_threshold).fillna(False)
    passed = (area_ok & total_ok).astype(bool)

    adata.obs["passed_QC"] = passed.values
    adata.obs["total_intensity"] = total.values
    if area_col is None:
        adata.obs["area"] = area.values
    adata.obs["qc_area_ok"] = area_ok.values.astype(bool)
    adata.obs["qc_total_ok"] = total_ok.values.astype(bool)
    adata.obs["qc_filter_reason"] = _failure_reason(area_ok, total_ok)

    adata.uns["qc"] = {
        "method": "adaptive_quantile_qc_filter",
        "description": "Dataset-specific lower quantile thresholds on area and total protein intensity.",
        "area_quantile": area_quantile,
        "total_quantile": total_quantile,
        "area_threshold": area_threshold,
        "total_threshold": total_threshold,
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
