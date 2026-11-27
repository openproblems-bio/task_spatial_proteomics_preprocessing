"""Pixel-level signal-versus-background metrics. Thin wrapper over `sprotmetrics`.

The only metric component that reads pixels rather than the cell table. The
expensive pass over the image is `sprotmetrics.cell_pixel_stats`, done once
here so the reported thresholds and object counts describe the same pass the
metrics were computed from.
"""

import argparse
import json

import anndata as ad
import numpy as np
import sprotmetrics as spm
from sprotmetrics._signal import qc_mask_from_table, signal_background_from_stats

QC_COLUMN = "passed_QC"


def _load_markers(path, n_channels):
    if not path:
        return [f"marker_{i}" for i in range(n_channels)]
    with open(path) as handle:
        payload = json.load(handle)
    markers = payload["markers"] if isinstance(payload, dict) else payload
    markers = [str(m) for m in markers]
    if len(markers) != n_channels:
        raise ValueError(f"got {len(markers)} marker names for {n_channels} channels")
    return markers


def run(par: dict):
    image = np.load(par["input_image"])
    if image.ndim == 2:
        image = image[None, :, :]
    labels = np.load(par["input_labels"])
    markers = _load_markers(par.get("input_markers"), image.shape[0])

    dilation_iterations = int(par.get("dilation_iterations", 3))
    norm_percentile = float(par.get("norm_percentile", 95.0))
    random_state = int(par.get("random_state", 0))

    stats = spm.cell_pixel_stats(
        image,
        labels,
        markers,
        dilation_iterations=dilation_iterations,
        random_state=random_state,
    )

    # QC lives in the cell table, not in the pixels, so it has to be carried
    # across by label id.
    qc_mask = None
    if par.get("input_qc"):
        qc_mask = qc_mask_from_table(
            ad.read_h5ad(par["input_qc"]), stats.means.shape[0], QC_COLUMN
        )
    if qc_mask is None:
        qc_mask = np.zeros(stats.means.shape[0], dtype=bool)

    scores = signal_background_from_stats(
        stats, norm_percentile=norm_percentile, qc_mask=qc_mask
    )

    out = ad.AnnData(
        uns={
            "metric_ids": np.array(list(scores), dtype=object),
            "metric_values": np.array(list(scores.values()), dtype=float),
            "metric_group": "signal_background",
            "markers": np.array(markers, dtype=object),
            "otsu_thresholds": stats.thresholds.astype(float),
            "dilation_iterations": dilation_iterations,
            "norm_percentile": norm_percentile,
            "n_objects": int(stats.present.sum()),
            "n_background_samples": int(stats.n_background_pseudocells),
            "qc_column": QC_COLUMN,
        }
    )
    out.write_h5ad(par["output"], compression="gzip")


## VIASH START
par = {
    "input_image": "corrected_image.npy",
    "input_labels": "segmentation_labels.npy",
    "input_markers": "corrected_markers.json",
    "input_qc": None,
    "dilation_iterations": 3,
    "norm_percentile": 95.0,
    "random_state": 0,
    "output": "metrics_signal_background.h5ad",
}
meta = {"name": "signal_background"}
## VIASH END


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input_image", required=True)
    ap.add_argument("--input_labels", required=True)
    ap.add_argument("--input_markers", default=None)
    ap.add_argument("--input_qc", default=None)
    ap.add_argument("--dilation_iterations", type=int, default=3)
    ap.add_argument("--norm_percentile", type=float, default=95.0)
    ap.add_argument("--random_state", type=int, default=0)
    ap.add_argument("--output", required=True)
    run(vars(ap.parse_args()))


if __name__ == "__main__":
    if "par" in globals():
        run(par)
    else:
        main()
