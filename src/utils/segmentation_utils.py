"""Helpers shared by the image-derived methods_proteomics_segmentation components.

Shipped to each component as a resource and imported via
`sys.path.append(meta["resources_dir"])`. Kept dependency-light (numpy, plus
scikit-image for relabelling) so that every segmentation image can carry it.
"""

import json

import numpy as np


def _is_unset(value) -> bool:
    return value is None or str(value).strip() in {"", "null", "None"}


def parse_marker_list(value):
    """Split a comma-separated marker string into a list of stripped names.

    Viash passes `multiple: false` strings verbatim, but accept a list too so
    the helpers can be called from tests and notebooks directly.
    """
    if _is_unset(value):
        return []
    if isinstance(value, (list, tuple)):
        items = value
    else:
        items = str(value).split(",")
    return [str(v).strip() for v in items if str(v).strip()]


def load_marker_names(path, n_channels: int):
    """Read the marker sidecar written by image correction.

    Accepts `{"markers": [...]}` (what the image-correction components emit) or
    a bare list. The number of names must match the image's channel count, as a
    mismatch would silently select the wrong channel.
    """
    with open(path) as f:
        blob = json.load(f)
    markers = blob["markers"] if isinstance(blob, dict) else blob
    markers = [str(m) for m in markers]
    if len(markers) != n_channels:
        raise ValueError(
            f"Marker sidecar {path} lists {len(markers)} markers, "
            f"but the image has {n_channels} channels."
        )
    return markers


def _rescale(channel: np.ndarray) -> np.ndarray:
    """Rescale a channel to [0, 1] between its 1st and 99.8th percentile."""
    finite = channel[np.isfinite(channel)]
    if finite.size == 0:
        return np.zeros_like(channel, dtype=np.float32)
    lo, hi = np.percentile(finite, [1.0, 99.8])
    if hi <= lo:
        return np.zeros_like(channel, dtype=np.float32)
    return np.clip((channel - lo) / (hi - lo), 0.0, 1.0).astype(np.float32)


def resolve_marker_indices(names, markers):
    """Map marker names to channel indices, case-insensitively."""
    lookup = {}
    for i, m in enumerate(markers):
        lookup.setdefault(m.lower(), i)
    missing = [n for n in names if n.lower() not in lookup]
    if missing:
        raise KeyError(f"Markers {missing} not found. Available markers: {markers}")
    return [lookup[n.lower()] for n in names]


def select_channel(img, markers, names, fallback_index=None):
    """Return one 2D float32 channel for segmentation.

    - `names` given: the named markers. A single marker is returned unchanged;
      several are each percentile-rescaled and max-projected, so that one bright
      stain does not drown the others.
    - no names, `fallback_index` given: `img[fallback_index]` unchanged.
    - neither: None, so callers can treat an optional channel as absent.

    Returns `(channel, description)` where description records what was used,
    for the metadata JSON.
    """
    names = parse_marker_list(names)
    if names:
        if markers is None:
            raise ValueError(
                f"Markers {names} were requested by name, but no --input_markers "
                "sidecar was given to resolve them to channels."
            )
        idx = resolve_marker_indices(names, markers)
        if len(idx) == 1:
            channel = img[idx[0]].astype(np.float32)
        else:
            channel = np.max(np.stack([_rescale(img[i].astype(np.float32)) for i in idx]), axis=0)
        return channel, {"markers": [markers[i] for i in idx], "channels": idx}

    if fallback_index is None:
        return None, None

    fallback_index = int(fallback_index)
    if fallback_index < 0 or fallback_index >= img.shape[0]:
        raise ValueError(
            f"nuclei_channel={fallback_index} out of range for image with {img.shape[0]} channels"
        )
    desc = {"channels": [fallback_index]}
    if markers is not None:
        desc["markers"] = [markers[fallback_index]]
    return img[fallback_index].astype(np.float32), desc


def load_image_and_channels(par, use_membrane: bool = False):
    """Load `--input_image` and pick the nuclear (and optionally membrane) channel.

    Returns `(image, nuclear, membrane, channel_info)`; `membrane` is None when
    `use_membrane` is False or no membrane markers were named.
    """
    image = np.load(par["input_image"])
    if image.ndim != 3:
        raise ValueError(f"Expected image in channel-y-x format, got shape {image.shape}")

    markers = None
    if not _is_unset(par.get("input_markers")):
        markers = load_marker_names(par["input_markers"], image.shape[0])

    nuclear, nuc_desc = select_channel(
        image, markers, par.get("nuclear_markers"), fallback_index=par.get("nuclei_channel", 0) or 0
    )
    if not np.isfinite(nuclear).any():
        raise ValueError("Nuclear channel contains no finite values.")

    membrane, mem_desc = (None, None)
    if use_membrane:
        membrane, mem_desc = select_channel(image, markers, par.get("membrane_markers"))

    info = {"nuclear": nuc_desc}
    if use_membrane:
        info["membrane"] = mem_desc
    return image, nuclear, membrane, info


def _accepts_max_size(func) -> bool:
    import inspect

    return "max_size" in inspect.signature(func).parameters


def remove_small_objects(ar, min_size: int, connectivity: int = 1):
    """`skimage.morphology.remove_small_objects` across scikit-image versions.

    scikit-image 0.26 deprecated `min_size` in favour of `max_size`, and 2.0
    drops `min_size`. The images install scikit-image unpinned, so pick the
    keyword at runtime. In 0.26 the deprecated `min_size=n` already behaves
    as `max_size=n` (checked: a 64 px object is removed by `min_size=64`), so
    passing the same number keeps results identical to what each version
    computes today, only without the deprecation warning.
    """
    from skimage import morphology

    if _accepts_max_size(morphology.remove_small_objects):
        return morphology.remove_small_objects(ar, max_size=int(min_size), connectivity=connectivity)
    return morphology.remove_small_objects(ar, min_size=int(min_size), connectivity=connectivity)


def remove_small_holes(ar, area_threshold: int, connectivity: int = 1):
    """`skimage.morphology.remove_small_holes` across scikit-image versions.

    Same deprecation as `remove_small_objects`: `area_threshold=n` became
    `max_size=n` in 0.26, with the same result.
    """
    from skimage import morphology

    if _accepts_max_size(morphology.remove_small_holes):
        return morphology.remove_small_holes(ar, max_size=int(area_threshold), connectivity=connectivity)
    return morphology.remove_small_holes(ar, area_threshold=int(area_threshold), connectivity=connectivity)


def finalize_labels(labels) -> np.ndarray:
    """Renumber a label image to a gapless 1..n int32 range.

    Downstream aggregation counts unique labels, and the component tests rely
    on `max(labels) == n_objects`.
    """
    from skimage.segmentation import relabel_sequential

    labels = np.asarray(labels)
    if labels.ndim != 2:
        raise ValueError(f"Expected a 2D label image, got shape {labels.shape}")
    if labels.size and labels.min() < 0:
        raise ValueError("Label image contains negative values.")
    labels, _, _ = relabel_sequential(labels.astype(np.int64, copy=False))
    if labels.max() > np.iinfo(np.int32).max:
        raise ValueError("Too many objects for an int32 label image.")
    return labels.astype(np.int32)


def write_outputs(par, labels, method: str, params: dict, channel_info=None):
    """Save the label image and, if requested, the metadata JSON."""
    np.save(par["output_labels"], labels)

    if not _is_unset(par.get("output_metadata")):
        metadata = {"method": method, **params}
        if channel_info is not None:
            metadata["input_channels"] = channel_info
        metadata["n_objects"] = int(labels.max())
        metadata["labels_shape"] = list(labels.shape)
        with open(par["output_metadata"], "w") as f:
            json.dump(metadata, f, indent=2, default=_json_default)


def _json_default(o):
    if isinstance(o, np.generic):
        return o.item()
    raise TypeError(f"Object of type {type(o).__name__} is not JSON serializable")
