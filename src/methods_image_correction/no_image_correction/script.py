import json
import numpy as np
import spatialdata as sd


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


def run(par: dict):
    input_zarr = par["input_zarr"]
    image_key = par.get("image_key")
    rescale_max = _as_bool(par.get("rescale_max", False))
    output_image = par["output_image"]
    output_markers = par["output_markers"]

    sdata = sd.read_zarr(input_zarr)

    if image_key in [None, "null", ""]:
        image_key = _pick_first_key(sdata.images, "images")

    img_obj = sdata.images[image_key]
    img_np = _ensure_cyx(_to_numpy(img_obj)).astype(np.float32)
    markers = _get_channel_names(img_obj, img_np.shape[0])

    corrected = np.nan_to_num(img_np, nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)
    corrected[corrected < 0] = 0

    if rescale_max:
        for c in range(corrected.shape[0]):
            mx = float(np.max(corrected[c]))
            if mx > 0:
                corrected[c] = corrected[c] / mx

    np.save(output_image, corrected)

    with open(output_markers, "w") as f:
        json.dump(
            {
                "method": "no_image_correction",
                "image_key": image_key,
                "markers": markers,
                "background_percentile": None,
                "rescale_max": rescale_max,
            },
            f,
            indent=2,
        )


## VIASH START
par = {
    "input_zarr": "input.zarr",
    "image_key": None,
    "rescale_max": False,
    "output_image": "corrected_image.npy",
    "output_markers": "markers.json",
}
meta = {"name": "no_image_correction"}
## VIASH END

if __name__ == "__main__":
    run(par)
