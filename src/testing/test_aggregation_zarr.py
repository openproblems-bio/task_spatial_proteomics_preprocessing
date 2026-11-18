"""Smoke test for the SpatialData-based protein aggregation component.

This is the older aggregation interface, which reads the image and the label
mask straight from a SpatialData store instead of from .npy arrays. The
benchmark workflow uses the array-based components (see test_aggregation.py);
this one is covered so that both interfaces stay runnable.
"""

import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np

## VIASH START
meta = {
    "name": "basic_protein_aggregation",
    "executable": "target/executable/methods_protein_aggregation/basic_protein_aggregation/basic_protein_aggregation",
    "resources_dir": "src/testing",
}
## VIASH END

sys.path.insert(0, meta["resources_dir"])
import isp_fixtures  # noqa: E402

name = meta["name"]
work = Path("test_work")
work.mkdir(exist_ok=True)

zarr_path = work / "fixture.zarr"
output_path = work / "cell_by_protein.h5ad"

print(">> Building a SpatialData fixture", flush=True)
isp_fixtures.make_spatialdata_zarr(zarr_path)
labels = isp_fixtures.build_labels()
n_objects = int(labels.max())

print(f">> Running {name}", flush=True)
subprocess.run(
    [
        meta["executable"],
        "--input_zarr", str(zarr_path),
        "--image_key", "fixture_image",
        "--labels_key", "fixture_cell",
        "--stat", "mean",
        "--output", str(output_path),
    ],
    check=True,
)

print(">> Checking output file exists", flush=True)
assert output_path.exists(), "cell_by_protein.h5ad was not written"

adata = ad.read_h5ad(output_path)

print(">> Checking the cell-by-protein contract", flush=True)
assert adata.shape == (n_objects, len(isp_fixtures.MARKERS)), (
    f"expected {n_objects} cells x {len(isp_fixtures.MARKERS)} markers, got {adata.shape}"
)
assert list(adata.var_names) == isp_fixtures.MARKERS, (
    f"marker names were not preserved: {list(adata.var_names)}"
)
assert "counts" in adata.layers, "layers['counts'] is required by the downstream components"
assert np.isfinite(adata.X).all(), "aggregated intensities contain NaN or infinite values"
assert (adata.X >= 0).all(), "aggregated intensities should not be negative"

print(">> Checking uns['aggregation']", flush=True)
aggregation = adata.uns["aggregation"]
assert aggregation["method"] == name, (
    f"uns['aggregation']['method'] is {aggregation['method']}, expected {name}"
)
assert aggregation["statistic"] == "mean"

print(">> Checking the aggregated values against the fixture", flush=True)
image = isp_fixtures.build_image()
expected = image[0][labels == 1].mean()
assert np.isclose(adata.X[0, 0], expected, rtol=1e-4), (
    f"mean intensity for cell 1, marker 0 is {adata.X[0, 0]}, expected {expected}"
)

print(">> All checks passed", flush=True)
