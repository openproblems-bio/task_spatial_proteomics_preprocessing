import numpy as np


def segment_binning(shape, bin_size: int) -> np.ndarray:
    """Square bins numbered row-major from 1.

    Equivalent to `txsim.preprocessing.segment_binning` (theislab/txsim, MIT),
    rewritten with integer arithmetic so the labels are exact integers.
    """
    if bin_size < 1:
        raise ValueError(f"bin_size must be a positive integer, got {bin_size}.")
    n_rows, n_cols = shape
    bins_per_row = -(-n_cols // bin_size)  # ceil
    row_bin = (np.arange(n_rows) // bin_size)[:, None]
    col_bin = (np.arange(n_cols) // bin_size)[None, :]
    return row_bin * bins_per_row + col_bin + 1


def run(par: dict):
    image = np.load(par["input_image"], mmap_mode="r")
    if image.ndim != 3:
        raise ValueError(f"Expected image in channel-y-x format, got shape {image.shape}")

    bin_size = int(par["bin_size"]) if par.get("bin_size") is not None else 30
    labels = finalize_labels(segment_binning(image.shape[1:], bin_size))

    write_outputs(par, labels, "binning_segmentation", {"bin_size": bin_size})


## VIASH START
par = {
    "input_image": "corrected_image.npy",
    "bin_size": 30,
    "output_labels": "segmentation_labels.npy",
    "output_metadata": None,
}
meta = {"name": "binning_segmentation", "resources_dir": "src/utils"}
## VIASH END

import sys  # noqa: E402

sys.path.append(meta["resources_dir"])
from segmentation_utils import finalize_labels, write_outputs  # noqa: E402

if __name__ == "__main__":
    run(par)
