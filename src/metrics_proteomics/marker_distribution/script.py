"""Per-marker distribution summaries.

Thin wrapper over `sprotmetrics.marker_distribution`; the per-marker table it
also publishes to `uns` comes from `sprotmetrics.marker_table`, which is the
same computation the summaries reduce.
"""

import argparse

import anndata as ad
import numpy as np
import sprotmetrics as spm

QC_COLUMN = "passed_QC"


def run(par: dict):
    layer = par.get("layer", "counts")
    zero_threshold = float(par.get("zero_threshold", 0.0))
    coverage_threshold = float(par.get("coverage_threshold", 0.1))
    low_q = float(par.get("low_quantile", 0.01))
    high_q = float(par.get("high_quantile", 0.99))

    adata = ad.read_h5ad(par["input"])
    shared = {
        "layer": layer,
        "zero_threshold": zero_threshold,
        "coverage_threshold": coverage_threshold,
        "quantiles": (low_q, high_q),
    }
    scores = spm.marker_distribution(adata, qc_variants=True, qc_column=QC_COLUMN, **shared)
    marker_table = spm.marker_table(adata, **shared)

    uns = {
        "metric_ids": np.array(list(scores), dtype=object),
        "metric_values": np.array(list(scores.values()), dtype=float),
        "metric_group": "marker_distribution",
        "evaluated_layer": layer,
        "zero_threshold": zero_threshold,
        "coverage_threshold": coverage_threshold,
        "low_quantile": low_q,
        "high_quantile": high_q,
        "qc_column": QC_COLUMN,
        "marker_table": marker_table.reset_index(drop=True).to_dict(orient="list"),
    }
    if QC_COLUMN in adata.obs.columns and np.asarray(adata.obs[QC_COLUMN]).any():
        qc_table = spm.marker_table(adata, subset=QC_COLUMN, **shared)
        uns["marker_table_qc_pass"] = qc_table.reset_index(drop=True).to_dict(orient="list")

    ad.AnnData(uns=uns).write_h5ad(par["output"], compression="gzip")


## VIASH START
par = {
    "input": "cell_by_protein_norm.h5ad",
    "layer": "counts",
    "zero_threshold": 0.0,
    "coverage_threshold": 0.1,
    "low_quantile": 0.01,
    "high_quantile": 0.99,
    "output": "metrics_marker_distribution.h5ad",
}
meta = {"name": "marker_distribution"}
## VIASH END


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--layer", default="counts")
    ap.add_argument("--zero_threshold", type=float, default=0.0)
    ap.add_argument("--coverage_threshold", type=float, default=0.1)
    ap.add_argument("--low_quantile", type=float, default=0.01)
    ap.add_argument("--high_quantile", type=float, default=0.99)
    ap.add_argument("--output", required=True)
    run(vars(ap.parse_args()))


if __name__ == "__main__":
    if "par" in globals():
        run(par)
    else:
        main()
