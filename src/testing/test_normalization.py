"""Smoke test shared by every methods_proteomics_normalization component."""

import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np

## VIASH START
meta = {
    "name": "total_log1p_normalization",
    "executable": "target/executable/methods_proteomics_normalization/total_log1p_normalization/total_log1p_normalization",
    "resources_dir": "src/testing",
}
## VIASH END

sys.path.insert(0, meta["resources_dir"])
import isp_fixtures  # noqa: E402

name = meta["name"]
work = Path("test_work")
work.mkdir(exist_ok=True)

input_path = work / "cell_by_protein_qc.h5ad"
output_path = work / "cell_by_protein_norm.h5ad"

SCALE = 10000.0
COFACTOR = 5.0

original = isp_fixtures.make_cell_by_protein(input_path)
raw = np.asarray(original.layers["counts"], dtype=np.float64)

print(f">> Running {name}", flush=True)
subprocess.run(
    [
        meta["executable"],
        "--input", str(input_path),
        "--scale", str(SCALE),
        "--cofactor", str(COFACTOR),
        "--output", str(output_path),
    ],
    check=True,
)

print(">> Checking output file exists", flush=True)
assert output_path.exists(), "cell_by_protein_norm.h5ad was not written"

adata = ad.read_h5ad(output_path)

print(">> Checking the normalization contract", flush=True)
assert adata.shape == original.shape, f"normalization must not change shape: {adata.shape}"
assert "counts" in adata.layers, "layers['counts'] must be preserved for normalization diagnostics"
assert "normalized" in adata.layers, "layers['normalized'] is required by the normalization API"
assert np.allclose(adata.layers["counts"], original.layers["counts"]), (
    "the raw counts layer must survive normalization unchanged"
)

norm = np.asarray(adata.layers["normalized"], dtype=np.float64)
assert norm.shape == raw.shape
assert np.isfinite(norm).all(), "normalized values contain NaN or infinite values"

print(">> Checking uns['normalization'] metadata", flush=True)
assert adata.uns["normalization"]["method"] == name, (
    f"uns['normalization']['method'] is {adata.uns['normalization']['method']}, expected {name}"
)

print(">> Checking the transformation itself", flush=True)
if name == "no_normalization":
    assert np.allclose(norm, raw, rtol=1e-5), (
        "the no-normalization control must copy raw intensities through unchanged"
    )
elif name == "total_log1p_normalization":
    totals = raw.sum(axis=1, keepdims=True)
    totals[totals <= 0] = 1.0
    expected = np.log1p(raw / totals * SCALE)
    assert np.allclose(norm, expected, rtol=1e-4, atol=1e-5), (
        "normalized layer does not match per-cell total scaling followed by log1p"
    )
elif name == "normalize_total_intensity":
    # Superseded combined component; defaults to the log1p transform.
    totals = raw.sum(axis=1, keepdims=True)
    totals[totals <= 0] = 1.0
    expected = np.log1p(raw / totals * SCALE)
    assert np.allclose(norm, expected, rtol=1e-4, atol=1e-5), (
        "normalized layer does not match per-cell total scaling followed by log1p"
    )
elif name == "arcsinh_normalization":
    expected = np.arcsinh(raw / COFACTOR)
    assert np.allclose(norm, expected, rtol=1e-4, atol=1e-5), (
        "normalized layer does not match arcsinh(x / cofactor)"
    )
elif name == "clr_normalization":
    clipped = np.clip(raw, 0.0, None)
    geo = np.exp(np.log1p(clipped).mean(axis=1, keepdims=True))
    expected = np.log1p(clipped / geo)
    assert np.allclose(norm, expected, rtol=1e-4, atol=1e-5), (
        "normalized layer does not match per-cell CLR log1p(x / exp(mean(log1p(x))))"
    )
else:
    raise AssertionError(f"no expected transformation defined for component '{name}'")

if name != "no_normalization":
    assert not np.allclose(norm, raw), f"{name} left the intensities unchanged"

print(">> All checks passed", flush=True)
