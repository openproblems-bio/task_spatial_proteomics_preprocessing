"""Smoke test shared by every methods_protein_aggregation component."""

import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np

## VIASH START
meta = {
    "name": "mean_intensity_aggregation",
    "executable": "target/executable/methods_protein_aggregation/mean_intensity_aggregation/mean_intensity_aggregation",
    "resources_dir": "src/testing",
}
## VIASH END

sys.path.insert(0, meta["resources_dir"])
import isp_fixtures  # noqa: E402

name = meta["name"]
work = Path("test_work")
work.mkdir(exist_ok=True)

image_path = work / "corrected_image.npy"
labels_path = work / "segmentation_labels.npy"
markers_path = work / "corrected_markers.json"
output_path = work / "cell_by_protein.h5ad"

image = isp_fixtures.make_image(image_path)
labels = isp_fixtures.make_labels(labels_path)
markers = isp_fixtures.make_markers_json(markers_path)
n_objects = int(labels.max())

print(f">> Running {name}", flush=True)
subprocess.run(
    [
        meta["executable"],
        "--input_image", str(image_path),
        "--input_labels", str(labels_path),
        "--input_markers", str(markers_path),
        "--output", str(output_path),
    ],
    check=True,
)

print(">> Checking output file exists", flush=True)
assert output_path.exists(), "cell_by_protein.h5ad was not written"

print(">> Checking the cell-by-protein contract", flush=True)
adata = ad.read_h5ad(output_path)
assert adata.shape == (n_objects, len(markers)), (
    f"expected {n_objects} cells x {len(markers)} markers, got {adata.shape}"
)
assert list(adata.var_names) == markers, f"marker names were not preserved: {list(adata.var_names)}"
assert "counts" in adata.layers, "layers['counts'] is required by the downstream QC and metric components"
assert np.allclose(adata.layers["counts"], adata.X), "layers['counts'] should mirror X on aggregation output"
assert np.isfinite(adata.X).all(), "aggregated intensities contain NaN or infinite values"
assert (adata.X >= 0).all(), "aggregated intensities should not be negative"

for col in ["cell_id", "area", "centroid_y", "centroid_x"]:
    assert col in adata.obs.columns, f"obs['{col}'] is required by downstream components"

print(">> Checking areas against the fixture", flush=True)
expected_areas = np.bincount(labels.ravel())[1:]
assert np.array_equal(adata.obs["area"].to_numpy().astype(int), expected_areas), (
    "reported cell areas do not match the fixture label mask"
)

print(">> Checking the aggregation statistic", flush=True)
statistic = adata.uns["aggregation"]["statistic"]
assert adata.uns["aggregation"]["method"] == name

# Recompute the expected value for the first object and first channel directly
# from the fixture, so a mislabelled statistic cannot pass.
pixels = image[0][labels == 1]
expected = {"mean": np.mean, "sum": np.sum, "median": np.median}[statistic](pixels)
assert np.isclose(adata.X[0, 0], expected, rtol=1e-4), (
    f"{name} reports statistic '{statistic}' but the value for cell 1, marker 0 "
    f"is {adata.X[0, 0]}, expected {expected}"
)

print(">> Rejecting an image/labels shape mismatch", flush=True)
bad_labels_path = work / "bad_labels.npy"
np.save(bad_labels_path, labels[:32, :32])
bad = subprocess.run(
    [
        meta["executable"],
        "--input_image", str(image_path),
        "--input_labels", str(bad_labels_path),
        "--input_markers", str(markers_path),
        "--output", str(work / "bad.h5ad"),
    ],
    capture_output=True,
    text=True,
)
assert bad.returncode != 0, "aggregating an image and label mask of different shapes must fail"

print(">> All checks passed", flush=True)
