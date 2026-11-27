import json

import numpy as np
from skimage.filters import threshold_otsu
from skimage.measure import label
from skimage.morphology import binary_closing, disk


def _nuclear_threshold(values: np.ndarray, method: str, quantile: float) -> float:
    """Threshold the nuclear channel using the requested strategy.

    Unknown strategies raise instead of silently falling back, so that a
    mistyped --threshold_method is not reported as a successful run.
    """
    method = str(method).lower()
    if method == "quantile":
        return float(np.quantile(values, quantile))
    if method == "otsu":
        return float(threshold_otsu(values))
    raise ValueError(f"threshold_method must be either 'quantile' or 'otsu', got '{method}'.")


def run(par: dict):
    nuclei_channel = int(par.get("nuclei_channel", 0))
    threshold_method = par.get("threshold_method") or "quantile"
    threshold_quantile = float(par.get("threshold_quantile", 0.90))
    min_area = int(par.get("min_area", 30))
    output_labels = par["output_labels"]
    output_metadata = par.get("output_metadata")

    # Selects the nuclear channel by marker name when --nuclear_markers is set,
    # and by --nuclei_channel otherwise.
    _, nuc, _, channel_info = load_image_and_channels(par)
    finite = nuc[np.isfinite(nuc)]

    thr = _nuclear_threshold(finite, threshold_method, threshold_quantile)

    mask = nuc > thr
    mask = binary_closing(mask, disk(1))
    mask = remove_small_objects(mask, min_size=min_area)

    # label() on a boolean mask always numbers objects 1..n without gaps.
    labels = label(mask).astype(np.int32)

    np.save(output_labels, labels)

    if output_metadata:
        n_objects = int(labels.max())
        with open(output_metadata, "w") as f:
            json.dump(
                {
                    "method": "basic_threshold_segmentation",
                    "nuclei_channel": nuclei_channel,
                    "input_channels": channel_info,
                    "threshold_method": str(threshold_method).lower(),
                    "threshold_quantile": threshold_quantile,
                    "threshold_value": thr,
                    "min_area": min_area,
                    "n_objects": n_objects,
                    "labels_shape": list(labels.shape),
                },
                f,
                indent=2,
            )


## VIASH START
par = {
    "input_image": "corrected_image.npy",
    "input_markers": None,
    "nuclear_markers": None,
    "nuclei_channel": 0,
    "threshold_method": "quantile",
    "threshold_quantile": 0.90,
    "min_area": 30,
    "output_labels": "segmentation_labels.npy",
    "output_metadata": None,
}
meta = {"name": "basic_threshold_segmentation", "resources_dir": "src/utils"}
## VIASH END

import sys  # noqa: E402

sys.path.append(meta["resources_dir"])
from segmentation_utils import load_image_and_channels, remove_small_objects  # noqa: E402

if __name__ == "__main__":
    run(par)
