"""Dataset and object statistics.

The metric itself lives in `sprotmetrics`, which this component now depends on;
what remains here is the Viash contract: read the inputs, call the library,
write `uns['metric_ids']` and `uns['metric_values']` in the established order.
"""

import argparse

import anndata as ad
import numpy as np
import sprotmetrics as spm

QC_COLUMN = "passed_QC"


def run(par: dict):
    cov_thr = par.get("cov_thr", par.get("coverage_threshold", 0.1))
    cov_thr = 0.1 if cov_thr is None else float(cov_thr)

    adata = ad.read_h5ad(par["input"])
    scores = spm.basic_statistics(
        adata,
        layer="counts",
        area="area",
        threshold=cov_thr,
        qc_variants=True,
        qc_column=QC_COLUMN,
    )

    out = ad.AnnData(
        uns={
            "metric_ids": np.array(list(scores), dtype=object),
            "metric_values": np.array(list(scores.values()), dtype=float),
            "metric_group": "basic",
            "coverage_threshold": cov_thr,
            "qc_column": QC_COLUMN,
        }
    )
    out.write_h5ad(par["output"], compression="gzip")


## VIASH START
par = {
    "input": "cell_by_protein_norm.h5ad",
    "cov_thr": 0.1,
    "output": "metrics_basic.h5ad",
}
meta = {"name": "basic"}
## VIASH END


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--cov_thr", type=float, default=None)
    ap.add_argument("--coverage_threshold", type=float, default=None)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()
    cov_thr = args.cov_thr if args.cov_thr is not None else args.coverage_threshold
    run({"input": args.input, "cov_thr": cov_thr, "output": args.output})


if __name__ == "__main__":
    if "par" in globals():
        run(par)
    else:
        main()
