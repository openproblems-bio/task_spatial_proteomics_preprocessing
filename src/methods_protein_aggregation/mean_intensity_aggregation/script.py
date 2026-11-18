import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from skimage.measure import regionprops_table

## VIASH START
par = {
    "input_image": "corrected_image.npy",
    "input_labels": "segmentation_labels.npy",
    "input_markers": "corrected_markers.json",
    "output": "cell_by_protein.h5ad",
}
## VIASH END


def _load_markers(path, n_channels):
    with open(path) as f:
        data = json.load(f)
    if isinstance(data, dict):
        markers = data.get("markers") or data.get("channels") or data.get("marker_names")
    else:
        markers = data
    if markers is None or len(markers) != n_channels:
        markers = [f"marker_{i}" for i in range(n_channels)]
    return [str(x) for x in markers]


def _ensure_cyx(image):
    image = np.asarray(image)
    if image.ndim != 3:
        raise ValueError(f"Expected 3D image array, got shape {image.shape}.")
    return image


def _obs_from_labels(labels):
    props = regionprops_table(
        labels.astype(np.int64),
        properties=("label", "area", "centroid"),
    )
    obs = pd.DataFrame(props)
    obs = obs.rename(columns={"centroid-0": "centroid_y", "centroid-1": "centroid_x"})
    obs["cell_id"] = obs["label"].astype(str)
    obs.index = obs["cell_id"].astype(str)
    return obs


def _aggregate_mean(image, labels, label_ids):
    flat_labels = labels.ravel().astype(np.int64)
    valid = flat_labels > 0
    labs = flat_labels[valid]
    areas = np.bincount(labs, minlength=int(flat_labels.max()) + 1).astype(float)

    values = []
    for c in range(image.shape[0]):
        vals = image[c].ravel()[valid].astype(float)
        sums = np.bincount(labs, weights=vals, minlength=len(areas))
        means = sums[label_ids] / np.maximum(areas[label_ids], 1.0)
        values.append(means)
    return np.vstack(values).T


def run(par):
    image = _ensure_cyx(np.load(par["input_image"]))
    labels = np.asarray(np.load(par["input_labels"])).astype(np.int64)

    if image.shape[1:] != labels.shape:
        raise ValueError(f"Image spatial shape {image.shape[1:]} does not match labels shape {labels.shape}.")

    markers = _load_markers(par["input_markers"], image.shape[0])
    obs = _obs_from_labels(labels)
    label_ids = obs["label"].astype(int).to_numpy()
    X = _aggregate_mean(image, labels, label_ids)

    adata = ad.AnnData(
        X=X.astype(np.float32),
        obs=obs.drop(columns=["label"]),
        var=pd.DataFrame(index=markers),
    )
    adata.layers["counts"] = adata.X.copy()
    adata.uns["aggregation"] = {
        "method": "mean_intensity_aggregation",
        "statistic": "mean",
        "input_image": str(par["input_image"]),
        "input_labels": str(par["input_labels"]),
        "n_labels": int(len(obs)),
        "n_markers": int(len(markers)),
    }

    Path(par["output"]).parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(par["output"])


if __name__ == "__main__":
    run(par)
