"""Smoke test shared by every methods_proteomics_qc component."""

import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np

## VIASH START
meta = {
    "name": "basic_qc_filter",
    "executable": "target/executable/methods_proteomics_qc/basic_qc_filter/basic_qc_filter",
    "resources_dir": "src/testing",
}
## VIASH END

sys.path.insert(0, meta["resources_dir"])
import isp_fixtures  # noqa: E402

name = meta["name"]
work = Path("test_work")
work.mkdir(exist_ok=True)

input_path = work / "cell_by_protein.h5ad"
output_path = work / "cell_by_protein_qc.h5ad"

original = isp_fixtures.make_cell_by_protein(input_path)

print(f">> Running {name}", flush=True)
subprocess.run(
    [
        meta["executable"],
        "--input", str(input_path),
        "--min_area", "50",
        "--min_total", "0.0",
        "--area_quantile", "0.05",
        "--total_quantile", "0.05",
        "--output", str(output_path),
    ],
    check=True,
)

print(">> Checking output file exists", flush=True)
assert output_path.exists(), "cell_by_protein_qc.h5ad was not written"

adata = ad.read_h5ad(output_path)

print(">> Checking QC annotates rather than filters", flush=True)
# The workflow contract is that QC keeps every cell and only adds annotations,
# so that downstream steps can choose to inspect or drop failing cells.
assert adata.n_obs == original.n_obs, (
    f"QC must not drop cells: {original.n_obs} in, {adata.n_obs} out"
)
assert adata.n_vars == original.n_vars
assert np.allclose(adata.layers["counts"], original.layers["counts"]), (
    "QC must not modify the counts layer"
)

print(">> Checking the QC annotation contract", flush=True)
for col in ["passed_QC", "total_intensity", "qc_area_ok", "qc_total_ok", "qc_filter_reason"]:
    assert col in adata.obs.columns, f"obs['{col}'] is required by the QC API"

assert adata.obs["passed_QC"].dtype == bool, "obs['passed_QC'] must be boolean"
expected_totals = np.asarray(original.layers["counts"]).sum(axis=1)
assert np.allclose(adata.obs["total_intensity"].to_numpy(), expected_totals, rtol=1e-4), (
    "obs['total_intensity'] does not match the summed counts layer"
)

print(">> Checking uns['qc'] metadata", flush=True)
qc = adata.uns["qc"]
assert qc["method"] == name, f"uns['qc']['method'] is {qc['method']}, expected {name}"
assert qc["n_cells"] == adata.n_obs
assert qc["n_passed"] == int(adata.obs["passed_QC"].sum())
assert np.isclose(qc["pass_rate"], float(adata.obs["passed_QC"].mean()))

print(">> Checking method-specific behaviour", flush=True)
passed = adata.obs["passed_QC"].to_numpy()
reasons = set(adata.obs["qc_filter_reason"].unique())

if name == "no_qc_filter":
    assert passed.all(), "the no-QC control must pass every cell"
    assert reasons == {"pass_no_qc"}, f"unexpected QC reasons for the control: {reasons}"
else:
    # The fixture spans areas 5..400, so a min_area of 50 and a 5% area
    # quantile both have to reject something. A method that passes everything
    # here is not applying its thresholds.
    assert not passed.all(), f"{name} passed every cell; thresholds are not being applied"
    assert passed.any(), f"{name} rejected every cell on the fixture"
    assert "pass" in reasons, "passing cells should be labelled 'pass'"
    # Every failing cell must carry a reason other than 'pass'.
    failing_reasons = set(adata.obs.loc[~passed, "qc_filter_reason"].unique())
    assert failing_reasons and "pass" not in failing_reasons, (
        f"failing cells must state a reason, got {failing_reasons}"
    )
    # passed_QC must be exactly the conjunction of the two component flags.
    conjunction = adata.obs["qc_area_ok"].to_numpy() & adata.obs["qc_total_ok"].to_numpy()
    assert np.array_equal(passed, conjunction), (
        "obs['passed_QC'] should equal qc_area_ok & qc_total_ok"
    )

print(">> All checks passed", flush=True)
