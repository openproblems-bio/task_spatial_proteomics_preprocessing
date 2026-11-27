import json

import numpy as np
import spatialdata as sd


def _to_numpy(x):
    if hasattr(x, "data"):
        x = x.data
    if hasattr(x, "compute"):
        x = x.compute()
    return np.asarray(x)


def _as_2d_labels(label_obj) -> np.ndarray:
    arr = _to_numpy(label_obj)
    dims = list(getattr(label_obj, "dims", []))

    if dims and len(dims) == arr.ndim and "y" in dims and "x" in dims:
        y_axis = dims.index("y")
        x_axis = dims.index("x")
        index = [0] * arr.ndim
        index[y_axis] = slice(None)
        index[x_axis] = slice(None)
        arr = arr[tuple(index)]
    else:
        arr = np.squeeze(arr)
        if arr.ndim == 3 and arr.shape[0] == 1:
            arr = arr[0]
        elif arr.ndim == 3 and arr.shape[-1] == 1:
            arr = arr[..., 0]

    arr = np.squeeze(arr)
    if arr.ndim != 2:
        raise ValueError(f"Expected a 2D label image after squeezing/indexing, got shape {arr.shape}.")

    return arr.astype(np.int32, copy=False)


def run(par: dict):
    input_zarr = par.get("input_zarr")
    labels_key = par.get("labels_key", None)
    output_labels = par["output_labels"]
    output_metadata = par.get("output_metadata")

    if input_zarr is None or str(input_zarr).strip() == "":
        raise ValueError("reference_segmentation requires --input_zarr.")

    sdata = sd.read_zarr(input_zarr)

    if not hasattr(sdata, "labels") or len(sdata.labels) == 0:
        raise ValueError("No labels elements found in the SpatialData object.")

    if labels_key is None or str(labels_key).strip() in {"", "null", "None"}:
        labels_key = list(sdata.labels.keys())[0]

    if labels_key not in sdata.labels:
        available = list(sdata.labels.keys())
        raise KeyError(f"labels_key='{labels_key}' not found. Available labels keys: {available}")

    labels = _as_2d_labels(sdata.labels[labels_key])

    if labels.min() < 0:
        raise ValueError("Reference labels contain negative label values, which are not supported.")

    labels = labels.astype(np.int32, copy=False)
    np.save(output_labels, labels)

    if output_metadata:
        # Reference label ids are deliberately preserved as-is (not renumbered),
        # so that predictions stay comparable to the dataset's own cell ids.
        unique_labels = np.unique(labels)
        n_objects = int((unique_labels > 0).sum())
        with open(output_metadata, "w") as f:
            json.dump(
                {
                    "method": "reference_segmentation",
                    "input_zarr": str(input_zarr),
                    "labels_key": str(labels_key),
                    "n_objects": n_objects,
                    "max_label": int(labels.max()),
                    "relabeled": False,
                    "labels_shape": list(labels.shape),
                },
                f,
                indent=2,
            )


## VIASH START
par = {
    "input_zarr": "data.zarr",
    "labels_key": None,
    "output_labels": "segmentation_labels.npy",
    "output_metadata": None,
}
meta = {"name": "reference_segmentation"}
## VIASH END

if __name__ == "__main__":
    run(par)
