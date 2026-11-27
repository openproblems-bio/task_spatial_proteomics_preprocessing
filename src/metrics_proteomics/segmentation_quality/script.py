"""Reference-free segmentation quality. Thin wrapper over `sprotmetrics`.

The five metrics come from two places in the library: panel-knowledge checks,
and `spillover_ratio`, which needs centroids instead. The split is an internal
detail of the package; the ids and their order here are unchanged.
"""

import argparse

import anndata as ad
import numpy as np
import sprotmetrics as spm

QC_COLUMN = "passed_QC"


def _clean(value):
    """Treat Viash's 'null' placeholder as absent."""
    if value is None:
        return None
    text = str(value).strip()
    return None if text.lower() in {"", "null", "none"} else text


def run(par: dict):
    layer = par.get("layer", "counts")
    n_neighbors = int(par.get("n_neighbors", 6))
    max_cells = int(par.get("max_cells", 5000))
    random_state = int(par.get("random_state", 0))
    nuclear_markers = _clean(par.get("nuclear_markers"))
    exclusive_pairs = _clean(par.get("exclusive_pairs"))

    adata = ad.read_h5ad(par["input"])
    scores = spm.segmentation_quality(
        adata,
        layer=layer,
        nuclear_markers=nuclear_markers,
        exclusive_pairs=exclusive_pairs,
        n_neighbors=n_neighbors,
        max_cells=max_cells,
        random_state=random_state,
        qc_variants=True,
        qc_column=QC_COLUMN,
    )

    markers = [str(m) for m in adata.var_names]
    resolved_nuclear = [markers[i] for i in spm.resolve_nuclear_markers(markers, nuclear_markers)]
    # Kept as "A:B" strings, which is the shape downstream readers already expect.
    resolved_pairs = [
        f"{markers[i]}:{markers[j]}"
        for i, j in spm.resolve_exclusive_pairs(markers, exclusive_pairs)
    ]

    out = ad.AnnData(
        uns={
            "metric_ids": np.array(list(scores), dtype=object),
            "metric_values": np.array(list(scores.values()), dtype=float),
            "metric_group": "segmentation_quality",
            "evaluated_layer": layer,
            "nuclear_markers": np.array(resolved_nuclear, dtype=object),
            "exclusive_pairs": np.array(resolved_pairs, dtype=object),
            "n_neighbors": n_neighbors,
            "qc_column": QC_COLUMN,
        }
    )
    out.write_h5ad(par["output"], compression="gzip")


## VIASH START
par = {
    "input": "cell_by_protein_norm.h5ad",
    "layer": "counts",
    "nuclear_markers": None,
    "exclusive_pairs": None,
    "n_neighbors": 6,
    "max_cells": 5000,
    "random_state": 0,
    "output": "metrics_segmentation_quality.h5ad",
}
meta = {"name": "segmentation_quality"}
## VIASH END


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input", required=True)
    ap.add_argument("--layer", default="counts")
    ap.add_argument("--nuclear_markers", default=None)
    ap.add_argument("--exclusive_pairs", default=None)
    ap.add_argument("--n_neighbors", type=int, default=6)
    ap.add_argument("--max_cells", type=int, default=5000)
    ap.add_argument("--random_state", type=int, default=0)
    ap.add_argument("--output", required=True)
    run(vars(ap.parse_args()))


if __name__ == "__main__":
    if "par" in globals():
        run(par)
    else:
        main()
