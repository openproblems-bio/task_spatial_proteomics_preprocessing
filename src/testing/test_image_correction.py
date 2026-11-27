"""Smoke test shared by every methods_image_correction component."""

import json
import subprocess
import sys
from pathlib import Path

import numpy as np

## VIASH START
meta = {
    "name": "basic_image_correction",
    "executable": "target/executable/methods_image_correction/basic_image_correction/basic_image_correction",
    "resources_dir": "src/testing",
}
## VIASH END

sys.path.insert(0, meta["resources_dir"])
import isp_fixtures  # noqa: E402

name = meta["name"]
work = Path("test_work")
work.mkdir(exist_ok=True)

zarr_path = work / "fixture.zarr"
image_path = work / "corrected_image.npy"
markers_path = work / "corrected_markers.json"

print(">> Building a SpatialData fixture", flush=True)
isp_fixtures.make_spatialdata_zarr(zarr_path)

print(f">> Running {name}", flush=True)
subprocess.run(
    [
        meta["executable"],
        "--input_zarr", str(zarr_path),
        "--image_key", "fixture_image",
        "--output_image", str(image_path),
        "--output_markers", str(markers_path),
    ],
    check=True,
)

print(">> Checking output files exist", flush=True)
assert image_path.exists(), "corrected_image.npy was not written"
assert markers_path.exists(), "corrected_markers.json was not written"

print(">> Checking corrected image", flush=True)
img = np.load(image_path)
assert img.shape == isp_fixtures.IMAGE_SHAPE, (
    f"expected channel-y-x shape {isp_fixtures.IMAGE_SHAPE}, got {img.shape}"
)
assert img.dtype == np.float32, f"expected float32, got {img.dtype}"
assert np.isfinite(img).all(), "corrected image contains NaN or infinite values"
assert img.min() >= 0.0, "corrected image contains negative intensities"
assert img.max() > 0.0, "corrected image is empty"

print(">> Checking marker metadata", flush=True)
with markers_path.open() as f:
    markers = json.load(f)
assert isinstance(markers, dict), "corrected_markers.json should be a JSON object"
assert markers["markers"] == isp_fixtures.MARKERS, (
    f"marker names were not preserved: {markers['markers']}"
)
assert markers["image_key"] == "fixture_image"
assert markers["method"] == name, f"metadata reports method {markers['method']}, expected {name}"

print(">> Checking --rescale_max", flush=True)
rescaled_path = work / "rescaled.npy"
subprocess.run(
    [
        meta["executable"],
        "--input_zarr", str(zarr_path),
        "--image_key", "fixture_image",
        "--rescale_max", "true",
        "--output_image", str(rescaled_path),
        "--output_markers", str(work / "rescaled_markers.json"),
    ],
    check=True,
)
rescaled = np.load(rescaled_path)
per_channel_max = rescaled.reshape(rescaled.shape[0], -1).max(axis=1)
assert np.allclose(per_channel_max, 1.0, atol=1e-5), (
    f"--rescale_max should leave every channel at max 1, got {per_channel_max}"
)

print(">> All checks passed", flush=True)
