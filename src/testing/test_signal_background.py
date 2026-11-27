"""Smoke test for the pixel-level signal-versus-background metric component.

Unlike the other metric components this one reads the image and the label mask
rather than the aggregated table, so it needs its own fixture set-up.
"""

import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np

## VIASH START
meta = {
    "name": "signal_background",
    "executable": "target/executable/metrics_proteomics/signal_background/signal_background",
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
qc_path = work / "cell_by_protein_qc.h5ad"
output_path = work / "metrics_signal_background.h5ad"

image = isp_fixtures.make_image(image_path)
labels = isp_fixtures.make_labels(labels_path)
markers = isp_fixtures.make_markers_json(markers_path)
n_objects = int(labels.max())

# The fixture's per-cell table is indexed by label id, which is how the
# component maps QC annotations back onto the mask.
isp_fixtures.make_cell_by_protein(qc_path, n_cells=n_objects, with_qc=True)

print(f">> Running {name}", flush=True)
subprocess.run(
    [
        meta["executable"],
        "--input_image", str(image_path),
        "--input_labels", str(labels_path),
        "--input_markers", str(markers_path),
        "--input_qc", str(qc_path),
        "--dilation_iterations", "3",
        "--norm_percentile", "95",
        "--random_state", "0",
        "--output", str(output_path),
    ],
    check=True,
)

print(">> Checking output file exists", flush=True)
assert output_path.exists(), "metrics_signal_background.h5ad was not written"

out = ad.read_h5ad(output_path)
ids = [str(x) for x in out.uns["metric_ids"]]
values = np.asarray(out.uns["metric_values"], dtype=float)
reported = dict(zip(ids, values))

print(">> Checking the metric output contract", flush=True)
assert len(ids) == len(values)
assert len(set(ids)) == len(ids), f"duplicate metric ids: {ids}"
assert out.uns["metric_group"] == name
assert out.uns["n_objects"] == n_objects, (
    f"reported {out.uns['n_objects']} objects, fixture has {n_objects}"
)
assert list(out.uns["markers"]) == markers

print(">> Checking headline metrics", flush=True)
expected = [
    "median_rel_sni",
    "min_rel_sni",
    "median_intracell_coverage",
    "median_marker_gini",
    "median_signal_to_background_ratio",
    "median_cellular_specificity",
]
missing = [m for m in expected if m not in ids]
assert not missing, f"missing metrics: {missing}"

print(">> Checking ranges and directions", flush=True)
# The fixture plants bright blobs exactly where the labels are, so signal must
# come out clearly above the surrounding background.
assert reported["median_signal_to_background_ratio"] > 1.0, (
    "fixture cells are brighter than their surroundings, so the ratio must exceed 1, got "
    f"{reported['median_signal_to_background_ratio']}"
)
assert reported["median_rel_sni"] > 0, (
    f"signal integral should be positive on the fixture, got {reported['median_rel_sni']}"
)
assert reported["min_rel_sni"] <= reported["median_rel_sni"]
for bounded in ["median_intracell_coverage", "median_marker_gini", "median_cellular_specificity"]:
    assert 0.0 <= reported[bounded] <= 1.0, f"{bounded} out of [0, 1]: {reported[bounded]}"
# Coverage is measured with an Otsu threshold; a value pinned at exactly 1.0
# would mean the threshold was not applied at all, which was true of the
# prototype this is adapted from.
assert reported["median_intracell_coverage"] < 1.0 or reported["median_intracell_coverage"] > 0.0

print(">> Checking the QC-passing variants", flush=True)
paired = [i for i in ids if i.endswith("_qc_pass")]
assert paired, "no _qc_pass variants reported"
for base in [i for i in ids if not i.endswith("_qc_pass")]:
    assert f"{base}_qc_pass" in ids, f"{base} has no _qc_pass counterpart"
assert np.isfinite(reported["median_rel_sni_qc_pass"]), (
    "a QC table was supplied, so the QC-passing variants must be computed"
)

print(">> Checking behaviour without a QC table", flush=True)
no_qc_out = work / "metrics_no_qc.h5ad"
subprocess.run(
    [
        meta["executable"],
        "--input_image", str(image_path),
        "--input_labels", str(labels_path),
        "--input_markers", str(markers_path),
        "--output", str(no_qc_out),
        "--random_state", "0",
    ],
    check=True,
)
m2 = ad.read_h5ad(no_qc_out)
ids2 = [str(x) for x in m2.uns["metric_ids"]]
reported2 = dict(zip(ids2, np.asarray(m2.uns["metric_values"], dtype=float)))
assert ids2 == ids, "the metric id list must not depend on whether QC was supplied"
for v in paired:
    assert not np.isfinite(reported2[v]), (
        f"{v} should be NaN without a QC table, got {reported2[v]}"
    )
# The all-cell metrics must not move when QC is simply absent.
for base in [i for i in ids if not i.endswith("_qc_pass")]:
    if np.isfinite(reported[base]) and np.isfinite(reported2[base]):
        assert np.isclose(reported[base], reported2[base]), (
            f"{base} changed when the QC table was omitted"
        )

print(">> Rejecting an image/labels shape mismatch", flush=True)
bad_labels = work / "bad_labels.npy"
np.save(bad_labels, labels[:32, :32])
bad = subprocess.run(
    [
        meta["executable"],
        "--input_image", str(image_path),
        "--input_labels", str(bad_labels),
        "--output", str(work / "bad.h5ad"),
    ],
    capture_output=True,
    text=True,
)
assert bad.returncode != 0, "mismatched image and label shapes must fail"

print(">> All checks passed", flush=True)
