import json
import numpy as np
import spatialdata as sd
from scipy.ndimage import gaussian_filter


def _pick_first_key(mapping, label: str):
    keys = list(mapping.keys())
    if not keys:
        raise ValueError(f"No {label} found in SpatialData.")
    return keys[0]


def _to_numpy(obj):
    return np.asarray(obj)


def _ensure_cyx(img_np: np.ndarray) -> np.ndarray:
    if img_np.ndim == 2:
        return img_np[None, :, :]

    if img_np.ndim == 3:
        s0, s1, s2 = img_np.shape
        if s0 <= 256 and s1 > 16 and s2 > 16:
            return img_np
        if s2 <= 256 and s0 > 16 and s1 > 16:
            return np.moveaxis(img_np, -1, 0)

    raise ValueError(f"Unsupported or ambiguous image shape: {img_np.shape}")


def _get_channel_names(image_obj, n_ch: int):
    for coord_name in ["c", "channel", "channels", "marker", "markers"]:
        coords = getattr(image_obj, "coords", None)
        if coords is not None and coord_name in coords:
            vals = coords[coord_name].values
            if len(vals) == n_ch:
                return [str(v) for v in vals]
    return [f"marker_{i}" for i in range(n_ch)]


def _as_bool(x) -> bool:
    if isinstance(x, bool):
        return x
    if x is None:
        return False
    return str(x).lower() in ["true", "1", "yes", "y"]


def _finite_percentile(x: np.ndarray, q: float) -> float:
    finite = x[np.isfinite(x)]
    if finite.size == 0:
        return 0.0
    return float(np.percentile(finite, q))


def _flatfield_correct_channel(
    x: np.ndarray,
    background_percentile: float,
    flatfield_sigma: float,
    flatfield_mode: str,
) -> np.ndarray:
    x = np.nan_to_num(x.astype(np.float32), nan=0.0, posinf=0.0, neginf=0.0)

    # Robust background subtraction before estimating the illumination field.
    bg = _finite_percentile(x, background_percentile)
    x0 = x - bg
    x0[x0 < 0] = 0

    # Smooth low-frequency field. This is a classical baseline, not a learned model.
    field = gaussian_filter(x0, sigma=float(flatfield_sigma))
    field = np.nan_to_num(field.astype(np.float32), nan=0.0, posinf=0.0, neginf=0.0)

    mode = str(flatfield_mode).lower()
    if mode == "subtract":
        # Remove low-frequency spatial bias while preserving the global intensity scale.
        baseline = float(np.median(field[field > 0])) if np.any(field > 0) else 0.0
        corrected = x0 - (field - baseline)
    elif mode == "divide":
        # Multiplicative flat-field correction. The median rescales values to the
        # original approximate intensity scale.
        positive = field[field > 0]
        baseline = float(np.median(positive)) if positive.size else 1.0
        eps = max(baseline * 1e-6, 1e-6)
        corrected = x0 / np.maximum(field, eps) * baseline
    else:
        raise ValueError("flatfield_mode must be either 'subtract' or 'divide'.")

    corrected = np.nan_to_num(corrected, nan=0.0, posinf=0.0, neginf=0.0)
    corrected[corrected < 0] = 0
    return corrected.astype(np.float32)


def run(par: dict):
    input_zarr = par["input_zarr"]
    image_key = par.get("image_key")
    background_percentile = float(par.get("background_percentile", 1.0))
    rescale_max = _as_bool(par.get("rescale_max", False))
    flatfield_sigma = float(par.get("flatfield_sigma", 50.0))
    flatfield_mode = str(par.get("flatfield_mode", "subtract"))
    output_image = par["output_image"]
    output_markers = par["output_markers"]

    sdata = sd.read_zarr(input_zarr)

    if image_key in [None, "null", ""]:
        image_key = _pick_first_key(sdata.images, "images")

    img_obj = sdata.images[image_key]
    img_np = _ensure_cyx(_to_numpy(img_obj)).astype(np.float32)
    markers = _get_channel_names(img_obj, img_np.shape[0])

    corrected = np.empty_like(img_np, dtype=np.float32)

    for c in range(img_np.shape[0]):
        x = _flatfield_correct_channel(
            img_np[c],
            background_percentile=background_percentile,
            flatfield_sigma=flatfield_sigma,
            flatfield_mode=flatfield_mode,
        )

        if rescale_max:
            mx = float(np.max(x))
            if mx > 0:
                x = x / mx

        corrected[c] = x

    np.save(output_image, corrected)

    with open(output_markers, "w") as f:
        json.dump(
            {
                "method": "flatfield_image_correction",
                "image_key": image_key,
                "markers": markers,
                "background_percentile": background_percentile,
                "flatfield_sigma": flatfield_sigma,
                "flatfield_mode": flatfield_mode,
                "rescale_max": rescale_max,
            },
            f,
            indent=2,
        )


## VIASH START
par = {
    "input_zarr": "input.zarr",
    "image_key": None,
    "background_percentile": 1.0,
    "flatfield_sigma": 50.0,
    "flatfield_mode": "subtract",
    "rescale_max": False,
    "output_image": "corrected_image.npy",
    "output_markers": "markers.json",
}
meta = {"name": "flatfield_image_correction"}
## VIASH END

if __name__ == "__main__":
    run(par)
