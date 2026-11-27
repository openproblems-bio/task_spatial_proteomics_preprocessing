import json

import numpy as np
from scipy import ndimage as ndi
from skimage.filters import gaussian, threshold_otsu
from skimage.feature import peak_local_max
from skimage.measure import label
from skimage.morphology import binary_closing, disk
from skimage.segmentation import watershed, expand_labels, relabel_sequential


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


def _make_markers(distance: np.ndarray, mask: np.ndarray, min_distance: int) -> np.ndarray:
    """Create watershed seed markers from distance-transform local maxima."""
    if mask.sum() == 0:
        return np.zeros(mask.shape, dtype=np.int32)

    coords = peak_local_max(
        distance,
        labels=mask,
        min_distance=max(1, int(min_distance)),
        exclude_border=False,
    )

    marker_img = np.zeros(mask.shape, dtype=bool)
    if coords.size > 0:
        marker_img[tuple(coords.T)] = True
        markers = label(marker_img).astype(np.int32)
    else:
        # Fallback: connected components if no robust local maxima are detected.
        markers = label(mask).astype(np.int32)

    return markers


def run(par: dict):
    nuclei_channel = int(par.get("nuclei_channel", 0))
    threshold_method = par.get("threshold_method") or "quantile"
    threshold_quantile = float(par.get("threshold_quantile", 0.90))
    min_area = int(par.get("min_area", 30))
    smoothing_sigma = float(par.get("smoothing_sigma", 1.0))
    hole_area = int(par.get("hole_area", 64))
    seed_min_distance = int(par.get("seed_min_distance", 5))
    expansion_distance = int(par.get("expansion_distance", 5))
    output_labels = par["output_labels"]
    output_metadata = par.get("output_metadata")

    # Selects the nuclear channel by marker name when --nuclear_markers is set,
    # and by --nuclei_channel otherwise.
    _, nuc, _, channel_info = load_image_and_channels(par)
    finite = nuc[np.isfinite(nuc)]

    if smoothing_sigma > 0:
        nuc_smooth = gaussian(nuc, sigma=smoothing_sigma, preserve_range=True)
    else:
        nuc_smooth = nuc

    thr = _nuclear_threshold(finite, threshold_method, threshold_quantile)
    mask = nuc_smooth > thr
    mask = binary_closing(mask, disk(1))
    mask = remove_small_holes(mask, area_threshold=max(1, hole_area))
    mask = remove_small_objects(mask, min_size=min_area)

    distance = ndi.distance_transform_edt(mask)
    markers = _make_markers(distance, mask, seed_min_distance)
    labels = watershed(-distance, markers=markers, mask=mask).astype(np.int32)
    labels = remove_small_objects(labels, min_size=min_area).astype(np.int32)

    if expansion_distance > 0 and labels.max() > 0:
        labels = expand_labels(labels, distance=expansion_distance).astype(np.int32)

    # Dropping small objects from a *label* image leaves holes in the numbering
    # (e.g. max label 7202 for only 6848 surviving objects). Downstream
    # aggregation counts unique labels, so renumber to a gapless 1..n range and
    # keep max(labels) == n_objects a usable invariant.
    labels, _, _ = relabel_sequential(labels)
    labels = labels.astype(np.int32)

    np.save(output_labels, labels)

    if output_metadata:
        n_objects = int(labels.max())
        with open(output_metadata, "w") as f:
            json.dump(
                {
                    "method": "watershed_expand_segmentation",
                    "nuclei_channel": nuclei_channel,
                    "input_channels": channel_info,
                    "threshold_method": str(threshold_method).lower(),
                    "threshold_quantile": threshold_quantile,
                    "threshold_value": thr,
                    "smoothing_sigma": smoothing_sigma,
                    "min_area": min_area,
                    "hole_area": hole_area,
                    "seed_min_distance": seed_min_distance,
                    "expansion_distance": expansion_distance,
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
    "smoothing_sigma": 1.0,
    "hole_area": 64,
    "seed_min_distance": 5,
    "expansion_distance": 5,
    "output_labels": "segmentation_labels.npy",
    "output_metadata": None,
}
meta = {"name": "watershed_expand_segmentation", "resources_dir": "src/utils"}
## VIASH END

import sys  # noqa: E402

sys.path.append(meta["resources_dir"])
from segmentation_utils import (  # noqa: E402
    load_image_and_channels,
    remove_small_holes,
    remove_small_objects,
)

if __name__ == "__main__":
    run(par)
