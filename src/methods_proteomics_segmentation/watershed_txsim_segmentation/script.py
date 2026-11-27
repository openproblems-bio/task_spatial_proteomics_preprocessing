"""Watershed segmentation ported from txsim.

The pipeline mirrors `txsim.preprocessing.segment_watershed` and its local
background filter (theislab/txsim, dev branch, MIT licence), as run by the
`watershed` component of openproblems-bio/task_ist_preprocessing. The generic
`*_func` dispatch and the None-valued pass-through arguments are reduced to the
choices txsim actually maps.
"""

import numpy as np
from scipy.ndimage import distance_transform_edt
from skimage import exposure, filters, morphology
from skimage.feature import peak_local_max
from skimage.filters import rank
from skimage.measure import label
from skimage.segmentation import relabel_sequential, watershed
from skimage.util import img_as_ubyte


def _is_none(value) -> bool:
    return value is None or str(value).strip().lower() in {"", "none", "null"}


def _footprint(name: str, size: int):
    name = str(name).lower()
    size = int(size)
    if name == "square":
        # skimage.morphology.square is deprecated; an all-ones block is identical.
        return np.ones((size, size), dtype=np.uint8)
    if name == "disk":
        return morphology.disk(size)
    if name == "diamond":
        return morphology.diamond(size)
    if name == "star":
        return morphology.star(size)
    raise ValueError(f"threshold_footprint must be square, disk, diamond or star, got '{name}'.")


def _normalize(img, par):
    func = str(par.get("normalize_func") or "none").lower()
    if func == "gamma":
        return exposure.adjust_gamma(img, gamma=par["normalize_gamma"], gain=par["normalize_gain"])
    if func == "log":
        return exposure.adjust_log(img, gain=par["normalize_gain"], inv=par["normalize_inv"])
    if func == "sigmoid":
        return exposure.adjust_sigmoid(
            img, cutoff=par["normalize_cutoff"], gain=par["normalize_gain"], inv=par["normalize_inv"]
        )
    if func == "none":
        return img
    raise ValueError(f"normalize_func must be gamma, log, sigmoid or none, got '{func}'.")


def _contrast(img, par):
    func = str(par.get("contrast_adjustment_func") or "none").lower()
    if func == "equalize_adapthist":
        kernel_size = par.get("contrast_adjustment_kernel_size")
        return exposure.equalize_adapthist(
            img,
            kernel_size=None if _is_none(kernel_size) else int(kernel_size),
            clip_limit=par["contrast_adjustment_clip_limit"],
            nbins=int(par["contrast_adjustment_nbins"]),
        )
    if func == "equalize_hist":
        return exposure.equalize_hist(img, nbins=int(par["contrast_adjustment_nbins"]))
    if func == "rescale_intensity":
        return exposure.rescale_intensity(img)
    if func == "none":
        return img
    raise ValueError(
        "contrast_adjustment_func must be equalize_adapthist, equalize_hist, "
        f"rescale_intensity or none, got '{func}'."
    )


def _blur(img, par):
    func = str(par.get("blur_func") or "none").lower()
    if func == "gaussian":
        return filters.gaussian(
            img,
            sigma=par["blur_sigma"],
            mode=par["blur_mode"],
            cval=par["blur_cval"],
            preserve_range=False,
            truncate=par["blur_truncate"],
        )
    if func == "median":
        return filters.median(img, mode=par["blur_mode"], cval=par["blur_cval"])
    if func == "none":
        return img
    raise ValueError(f"blur_func must be gaussian, median or none, got '{func}'.")


def _threshold(img_u8, par):
    func = str(par.get("threshold_func") or "otsu").lower()
    if func == "local_otsu":
        footprint = _footprint(par["threshold_footprint"], par["threshold_footprint_size"])
        return img_u8 >= rank.otsu(img_u8, footprint)
    if func == "otsu":
        return img_u8 > filters.threshold_otsu(img_u8, nbins=int(par["threshold_nbins"]))
    if func == "triangle":
        return img_u8 > filters.threshold_triangle(img_u8, nbins=int(par["threshold_nbins"]))
    raise ValueError(f"threshold_func must be local_otsu, otsu or triangle, got '{func}'.")


def _mean_intensity_per_cell(img, labels):
    """txsim get_mean_intensity_per_cell; labels must be gapless 1..n."""
    n = int(labels.max())
    sums = np.bincount(labels.ravel(), weights=img.ravel(), minlength=n + 1)[1:]
    sizes = np.bincount(labels.ravel(), minlength=n + 1)[1:].astype(float)
    out = np.zeros(n, dtype=float)
    np.divide(sums, sizes, out=out, where=sizes > 0)
    return out


def _local_background_per_cell(img, labels, window_size, bg_size):
    """txsim get_local_background_mean_intensity_per_cell.

    Each cell's background is the mean of unlabelled pixels in a `bg_size`
    window around the `window_size` tile containing it; a cell spanning tiles
    keeps the minimum.
    """
    if not window_size < bg_size:
        raise ValueError("bg_intensity_filter_window_size must be smaller than bg_intensity_filter_bg_size.")
    n = int(labels.max())
    bg = np.full(n + 1, float(np.max(img)))
    extend = int((bg_size - window_size) / 2)
    ny = -(-img.shape[0] // window_size)
    nx = -(-img.shape[1] // window_size)
    for i in range(nx):
        for j in range(ny):
            x0, x1 = i * window_size, min((i + 1) * window_size, img.shape[1])
            y0, y1 = j * window_size, min((j + 1) * window_size, img.shape[0])
            ids = np.unique(labels[y0:y1, x0:x1])
            ids = ids[ids > 0]
            if ids.size == 0:
                continue
            by0, by1 = max(0, y0 - extend), min(img.shape[0], y1 + extend)
            bx0, bx1 = max(0, x0 - extend), min(img.shape[1], x1 + extend)
            bg_pixels = img[by0:by1, bx0:bx1][labels[by0:by1, bx0:bx1] == 0]
            # As in txsim: a window with no background has a NaN background,
            # which propagates through np.minimum, and a NaN comparison never
            # filters the cell, so such cells are always kept.
            local = float(np.mean(bg_pixels)) if bg_pixels.size else np.nan
            bg[ids] = np.minimum(bg[ids], local)
    return bg[1:]


def _filter_by_local_background(img, labels, window_size, bg_size, bg_factor):
    """txsim filter_cells_based_on_local_background_intensity, with reindexing."""
    intensities = _mean_intensity_per_cell(img, labels)
    background = _local_background_per_cell(img, labels, window_size, bg_size)
    drop = np.concatenate([[False], intensities < bg_factor * background])
    labels = labels.copy()
    labels[drop[labels]] = 0
    labels, _, _ = relabel_sequential(labels)
    return labels


def _as_txsim_input(nuc: np.ndarray) -> np.ndarray:
    """Bring the nuclear channel into a range skimage's exposure functions accept.

    The IST images were unsigned integers, which txsim's pipeline consumes
    directly. Uncorrected proteomics images are integer-valued too (only
    stored as float), so those are cast back to uint16 losslessly and take
    exactly the txsim code path. Genuinely float images (after background
    subtraction or flat-field division) are clipped at zero and scaled to
    [0, 1], since img_as_uint/img_as_ubyte reject floats outside [-1, 1].
    """
    nuc = np.nan_to_num(nuc.astype(np.float64), nan=0.0, posinf=0.0, neginf=0.0)
    nuc = np.clip(nuc, 0.0, None)
    if nuc.max() <= np.iinfo(np.uint16).max and np.array_equal(nuc, np.round(nuc)):
        return nuc.astype(np.uint16)
    peak = nuc.max()
    return nuc / peak if peak > 0 else nuc


def segment_watershed(nuc: np.ndarray, par: dict) -> np.ndarray:
    nuc = _as_txsim_input(nuc)
    original = nuc.astype(np.float64)

    img = _normalize(nuc, par)
    img = _contrast(img, par)
    img = _blur(img, par)
    if np.issubdtype(img.dtype, np.floating) and (img.max() > 1 or img.min() < 0):
        # Only reachable on the float path with every step set to "none" or a
        # range-preserving option; integer images go to img_as_ubyte as in txsim.
        img = exposure.rescale_intensity(img, out_range=(0.0, 1.0))
    img_u8 = img_as_ubyte(img)

    foreground = _threshold(img_u8, par)
    connectivity = int(par["connectivity"])
    foreground = remove_small_objects(foreground, int(par["min_size"]), connectivity)
    foreground = remove_small_holes(foreground, int(par["hole_area_threshold"]), connectivity)

    distance = distance_transform_edt(foreground)
    peaks = peak_local_max(distance, min_distance=int(par["local_maxima_min_distance"]))
    marker_img = np.zeros(distance.shape, dtype=bool)
    if peaks.size:
        marker_img[tuple(peaks.T)] = True
    markers = label(marker_img)

    labels = watershed(
        -distance,
        markers,
        mask=foreground,
        watershed_line=True,
        compactness=float(par["watershed_compactness"]),
    )
    # txsim's background filter assumes gapless labels.
    labels, _, _ = relabel_sequential(labels)

    bg_factor = float(par["bg_intensity_filter_bg_factor"])
    if bg_factor > 0 and labels.max() > 0:
        labels = _filter_by_local_background(
            original,
            labels,
            window_size=int(par["bg_intensity_filter_window_size"]),
            bg_size=int(par["bg_intensity_filter_bg_size"]),
            bg_factor=bg_factor,
        )
    return labels


DEFAULTS = {
    "normalize_func": "gamma",
    "normalize_gamma": 1.0,
    "normalize_gain": 1.0,
    "normalize_inv": False,
    "normalize_cutoff": 0.5,
    "contrast_adjustment_func": "equalize_adapthist",
    "contrast_adjustment_kernel_size": None,
    "contrast_adjustment_clip_limit": 0.01,
    "contrast_adjustment_nbins": 256,
    "blur_func": "gaussian",
    "blur_sigma": 1.0,
    "blur_mode": "nearest",
    "blur_cval": 0.0,
    "blur_truncate": 4.0,
    "threshold_func": "local_otsu",
    "threshold_nbins": 256,
    "threshold_footprint": "square",
    "threshold_footprint_size": 50,
    "min_size": 64,
    "hole_area_threshold": 64,
    "connectivity": 1,
    "local_maxima_min_distance": 5,
    "watershed_compactness": 0.0,
    "bg_intensity_filter_bg_factor": 0.3,
    "bg_intensity_filter_window_size": 1000,
    "bg_intensity_filter_bg_size": 2000,
}


def run(par: dict):
    params = {k: (par.get(k) if par.get(k) is not None else v) for k, v in DEFAULTS.items()}

    _, nuc, _, channel_info = load_image_and_channels(par)
    labels = finalize_labels(segment_watershed(nuc, params))

    write_outputs(par, labels, "watershed_txsim_segmentation", params, channel_info)


## VIASH START
par = {
    "input_image": "corrected_image.npy",
    "input_markers": None,
    "nuclear_markers": None,
    "nuclei_channel": 0,
    "output_labels": "segmentation_labels.npy",
    "output_metadata": None,
}
meta = {"name": "watershed_txsim_segmentation", "resources_dir": "src/utils"}
## VIASH END

import sys  # noqa: E402

sys.path.append(meta["resources_dir"])
from segmentation_utils import (  # noqa: E402
    finalize_labels,
    load_image_and_channels,
    remove_small_holes,
    remove_small_objects,
    write_outputs,
)

if __name__ == "__main__":
    run(par)
