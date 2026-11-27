"""Smoke test shared by every methods_proteomics_segmentation component."""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np

## VIASH START
meta = {
    "name": "basic_threshold_segmentation",
    "executable": "target/executable/methods_proteomics_segmentation/basic_threshold_segmentation/basic_threshold_segmentation",
    "resources_dir": "src/testing",
}
## VIASH END

sys.path.insert(0, meta["resources_dir"])
import isp_fixtures  # noqa: E402

name = meta["name"]
work = Path("test_work")
work.mkdir(exist_ok=True)

labels_path = work / "segmentation_labels.npy"
metadata_path = work / "segmentation_metadata.json"

cmd = [meta["executable"], "--output_labels", str(labels_path), "--output_metadata", str(metadata_path)]

if name == "reference_segmentation":
    zarr_path = work / "fixture.zarr"
    isp_fixtures.make_spatialdata_zarr(zarr_path)
    cmd += ["--input_zarr", str(zarr_path), "--labels_key", "fixture_cell"]
else:
    image_path = work / "corrected_image.npy"
    isp_fixtures.make_image(image_path)
    # min_area is set above the smallest fixture blob so that object removal,
    # and therefore the label-renumbering path, is actually exercised.
    cmd += ["--input_image", str(image_path), "--nuclei_channel", "0", "--min_area", "20"]

print(f">> Running {name}", flush=True)
subprocess.run(cmd, check=True)

print(">> Checking output files exist", flush=True)
assert labels_path.exists(), "segmentation_labels.npy was not written"
assert metadata_path.exists(), "--output_metadata was requested but no file was written"

print(">> Checking label mask", flush=True)
labels = np.load(labels_path)
assert labels.ndim == 2, f"expected a 2D label mask, got shape {labels.shape}"
assert np.issubdtype(labels.dtype, np.integer), f"expected integer labels, got {labels.dtype}"
assert labels.min() >= 0, "label masks must not contain negative values"

present = np.unique(labels)
present = present[present > 0]
assert present.size > 0, "segmentation produced no objects on the fixture"

print(">> Checking metadata", flush=True)
with metadata_path.open() as f:
    metadata = json.load(f)
assert metadata["method"] == name, f"metadata reports method {metadata['method']}, expected {name}"
assert metadata["n_objects"] == present.size, (
    f"metadata n_objects={metadata['n_objects']} disagrees with "
    f"{present.size} distinct labels in the mask"
)
assert metadata["labels_shape"] == list(labels.shape)

if name == "reference_segmentation":
    # Reference ids are passed through unchanged so they stay comparable to the
    # dataset's own cell ids, so gapless numbering is not required here.
    assert metadata["relabeled"] is False
else:
    assert np.array_equal(present, np.arange(1, present.size + 1)), (
        f"labels must be a gapless 1..n range; got max={labels.max()} "
        f"for {present.size} objects"
    )
    assert int(labels.max()) == present.size, (
        f"max label ({labels.max()}) must equal the object count ({present.size})"
    )

print(f">> Rejecting an invalid --threshold_method", flush=True)
if name != "reference_segmentation":
    bad = subprocess.run(
        cmd + ["--threshold_method", "not_a_method"],
        capture_output=True,
        text=True,
    )
    assert bad.returncode != 0, "an unknown --threshold_method must fail instead of silently falling back"

print(">> All checks passed", flush=True)
