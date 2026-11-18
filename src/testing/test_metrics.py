"""Smoke test shared by every metrics_proteomics component."""

import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np

## VIASH START
meta = {
    "name": "basic",
    "executable": "target/executable/metrics_proteomics/basic/basic",
    "resources_dir": "src/testing",
}
## VIASH END

sys.path.insert(0, meta["resources_dir"])
import isp_fixtures  # noqa: E402

name = meta["name"]
work = Path("test_work")
work.mkdir(exist_ok=True)

input_path = work / "cell_by_protein_norm.h5ad"
output_path = work / f"metrics_{name}.h5ad"

# Metrics run at the end of the workflow, so the fixture carries everything the
# preceding steps would have added: counts, a normalized layer and QC columns.
# Both fixtures below must differ *only* in QC, or the comparison at the end of
# this script is meaningless. Keep this shared.
WITH_CELLTYPES = name in {"cluster_separability", "segmentation_quality"}
FIXTURE_CELLS = 200

original = isp_fixtures.make_cell_by_protein(
    input_path,
    # Large enough that the neighbour-based spillover metric has both adjacent
    # and non-adjacent negative cells to compare; at 40 cells every cell has a
    # positive neighbour and the metric is undefined.
    n_cells=FIXTURE_CELLS,
    with_normalized=True,
    with_qc=True,
    # cluster_separability scores obs["cell_type"]; segmentation_quality needs
    # the mutually exclusive lineages that the cell types create.
    with_celltypes=WITH_CELLTYPES,
)

# Metric components take different arguments; only --input and --output are shared.
extra_args = {
    "basic": ["--cov_thr", "0.1"],
    "marker_distribution": [
        "--layer", "counts",
        "--zero_threshold", "0.0",
        "--coverage_threshold", "0.1",
        "--low_quantile", "0.01",
        "--high_quantile", "0.99",
    ],
    "normalization_diagnostics": [
        "--raw_layer", "counts",
        "--normalized_layer", "normalized",
        "--max_cells", "5000",
        "--random_state", "0",
    ],
    "cluster_separability": [
        "--layer", "normalized",
        "--max_cells", "5000",
        "--random_state", "0",
    ],
    "segmentation_quality": [
        "--layer", "counts",
        "--nuclear_markers", ",".join(isp_fixtures.NUCLEAR_MARKERS),
        "--exclusive_pairs", ":".join(isp_fixtures.EXCLUSIVE_PAIR),
        "--n_neighbors", "4",
        "--random_state", "0",
    ],
    "technical_covariates": [
        "--layer", "normalized",
        "--covariate", "area",
        "--n_pcs", "3",
        "--random_state", "0",
    ],
}
assert name in extra_args, f"no argument set defined for metric component '{name}'"

print(f">> Running {name}", flush=True)
subprocess.run(
    [meta["executable"], "--input", str(input_path), "--output", str(output_path)] + extra_args[name],
    check=True,
)

print(">> Checking output file exists", flush=True)
assert output_path.exists(), f"metrics_{name}.h5ad was not written"

out = ad.read_h5ad(output_path)

print(">> Checking the metric output contract", flush=True)
assert "metric_ids" in out.uns, "uns['metric_ids'] is required by the metric API"
assert "metric_values" in out.uns, "uns['metric_values'] is required by the metric API"

ids = [str(x) for x in out.uns["metric_ids"]]
values = np.asarray(out.uns["metric_values"], dtype=float)

assert len(ids) == len(values), (
    f"uns['metric_ids'] has {len(ids)} entries but uns['metric_values'] has {len(values)}"
)
assert len(ids) > 0, "no metrics were reported"
assert len(set(ids)) == len(ids), f"duplicate metric ids reported: {ids}"
assert out.uns["metric_group"] == name, (
    f"uns['metric_group'] is {out.uns['metric_group']}, expected {name}"
)

print(">> Checking metric values are usable", flush=True)
# The collection script writes these straight into the summary table, so a NaN
# here silently becomes an empty cell in the published results.
nan_metrics = [i for i, v in zip(ids, values) if not np.isfinite(v)]
assert not nan_metrics, f"metrics are not finite on the fixture: {nan_metrics}"

print(">> Checking the family reports its headline metrics", flush=True)
expected_ids = {
    "basic": ["n_cells", "n_markers", "qc_pass_rate", "median_area", "median_total_intensity"],
    "marker_distribution": [
        "median_marker_cv",
        "median_marker_zero_fraction",
        "median_marker_coverage",
    ],
    "normalization_diagnostics": ["raw_total_cv", "normalized_total_cv", "total_cv_ratio_norm_over_raw"],
    "cluster_separability": ["silhouette_cell_type"],
    "segmentation_quality": [
        "dapi_positive_ratio",
        "replicate_correlation",
        "mutually_exclusive_coexpression_ratio",
        "spillover_ratio",
    ],
    "technical_covariates": [
        "variance_explained_by_covariate_normalized",
        "variance_explained_by_covariate_counts",
        "covariate_pc_correlation_normalized",
        "covariate_pc_correlation_counts",
        "mean_variance_correlation",
        "marker_variance_evenness",
    ],
}
missing = [m for m in expected_ids[name] if m not in ids]
assert not missing, f"{name} did not report expected metrics: {missing}"

reported = dict(zip(ids, values))

if name == "basic":
    assert reported["n_cells"] == original.n_obs, (
        f"n_cells is {reported['n_cells']}, expected {original.n_obs}"
    )
    assert reported["n_markers"] == original.n_vars
    expected_pass_rate = float(original.obs["passed_QC"].mean())
    assert np.isclose(reported["qc_pass_rate"], expected_pass_rate), (
        f"qc_pass_rate is {reported['qc_pass_rate']}, expected {expected_pass_rate}"
    )
    assert 0.0 <= reported["qc_pass_rate"] <= 1.0

if name == "marker_distribution":
    for bounded in ["median_marker_zero_fraction", "median_marker_coverage"]:
        assert 0.0 <= reported[bounded] <= 1.0, f"{bounded} is out of [0, 1]: {reported[bounded]}"

if name == "normalization_diagnostics":
    assert reported["raw_total_cv"] > 0, "the fixture has variable totals, so raw CV must be positive"
    # The fixture's normalized layer is total-scaled, so normalization must
    # compress the spread of per-cell totals rather than inflate it.
    assert reported["total_cv_ratio_norm_over_raw"] < 1.0, (
        "total-intensity CV should shrink after total-count normalization, got ratio "
        f"{reported['total_cv_ratio_norm_over_raw']}"
    )

if name == "cluster_separability":
    score = reported["silhouette_cell_type"]
    assert -1.0 <= score <= 1.0, f"silhouette must lie in [-1, 1], got {score}"
    # The fixture's two cell types differ in CD45/ECAD composition. The panel
    # also carries two correlated nuclear replicates, shared between the types,
    # which caps how separable they can look - so require clearly positive
    # rather than near-perfect.
    assert score > 0.3, f"expected clearly separated fixture clusters, got silhouette {score}"

if name == "segmentation_quality":
    # The fixture's two nuclear channels are the same signal restained, so
    # concordance must be high and per-cell disagreement small.
    assert reported["replicate_correlation"] > 0.8, (
        f"replicate_correlation should be high on replicate stains, got "
        f"{reported['replicate_correlation']}"
    )
    assert 0.0 <= reported["median_replicate_cv"] < 0.5
    assert 0.0 <= reported["dapi_positive_ratio"] <= 1.0
    # CD45-high and ECAD-high cells are disjoint by construction, so few cells
    # should be co-positive.
    assert 0.0 <= reported["mutually_exclusive_coexpression_ratio"] <= 1.0
    assert reported["mutually_exclusive_coexpression_ratio"] < 0.5, (
        "the fixture's lineages are disjoint, so co-expression should be low, got "
        f"{reported['mutually_exclusive_coexpression_ratio']}"
    )
    assert reported["spillover_ratio"] > 0

if name == "technical_covariates":
    bounded = ["marker_variance_evenness"]
    for suffix in ("normalized", "counts"):
        bounded += [
            f"variance_explained_by_covariate_{suffix}",
            f"max_marker_variance_explained_by_covariate_{suffix}",
            f"covariate_pc_correlation_{suffix}",
        ]
    for b in bounded:
        assert 0.0 <= reported[b] <= 1.0, f"{b} out of [0, 1]: {reported[b]}"
    assert -1.0 <= reported["mean_variance_correlation"] <= 1.0
    for suffix in ("normalized", "counts"):
        assert reported[f"max_marker_variance_explained_by_covariate_{suffix}"] >= reported[
            f"variance_explained_by_covariate_{suffix}"
        ], "the worst marker cannot explain less variance than the panel mean"

print(">> Checking the QC-passing variants", flush=True)
# QC components annotate rather than filter, so without a second pass over the
# QC-passing cells every metric here would be identical whichever QC method ran.
paired = [i for i in ids if i.endswith("_qc_pass")]
assert paired, f"{name} reported no _qc_pass variants"

# Metrics that cannot differ between the two cell sets carry no variant.
INVARIANT = {"n_markers", "qc_pass_rate"}
for base in ids:
    if base.endswith("_qc_pass") or base in INVARIANT:
        continue
    assert f"{base}_qc_pass" in ids, f"{base} has no _qc_pass counterpart"
for variant in paired:
    assert variant[: -len("_qc_pass")] in ids, f"{variant} has no base counterpart"

qc_pass = original.obs["passed_QC"].to_numpy()
assert not qc_pass.all(), "fixture must reject some cells for this check to mean anything"

if name == "basic":
    assert reported["n_cells_qc_pass"] == int(qc_pass.sum()), (
        f"n_cells_qc_pass is {reported['n_cells_qc_pass']}, "
        f"expected {int(qc_pass.sum())} passing cells"
    )
    # The fixture's rejected cells are debris: small and dim. Both summaries
    # must therefore rise once they are excluded.
    assert reported["median_area_qc_pass"] > reported["median_area"], (
        "median_area_qc_pass should exceed median_area when QC rejects small cells"
    )
    assert reported["median_total_intensity_qc_pass"] > reported["median_total_intensity"], (
        "median_total_intensity_qc_pass should exceed the all-cell value when QC "
        "rejects dim cells"
    )

# The two passes must actually disagree somewhere, or the subset is not being
# applied and the QC axis is still uninformative.
differing = [
    v for v in paired
    if np.isfinite(reported[v]) and np.isfinite(reported[v[: -len("_qc_pass")]])
    and not np.isclose(reported[v], reported[v[: -len("_qc_pass")]])
]
assert differing, (
    f"{name}: every _qc_pass metric equals its base value, so the QC subset "
    f"is not being applied"
)

print(">> Checking behaviour when QC has not been run", flush=True)
no_qc_path = work / "no_qc.h5ad"
isp_fixtures.make_cell_by_protein(
    no_qc_path,
    n_cells=FIXTURE_CELLS,
    with_normalized=True,
    with_qc=False,
    with_celltypes=WITH_CELLTYPES,
)
no_qc_out = work / f"metrics_{name}_no_qc.h5ad"
subprocess.run(
    [meta["executable"], "--input", str(no_qc_path), "--output", str(no_qc_out)] + extra_args[name],
    check=True,
)
m2 = ad.read_h5ad(no_qc_out)
ids2 = [str(x) for x in m2.uns["metric_ids"]]
vals2 = np.asarray(m2.uns["metric_values"], dtype=float)
reported2 = dict(zip(ids2, vals2))
assert ids2 == ids, "the metric id list must not depend on whether QC was run"
for v in paired:
    assert not np.isfinite(reported2[v]), (
        f"{v} should be NaN when obs['passed_QC'] is absent, got {reported2[v]}"
    )
# The all-cell metrics must be unaffected by the absence of QC annotations.
for base in ids:
    if base.endswith("_qc_pass") or base == "qc_pass_rate":
        continue
    if np.isfinite(reported[base]) and np.isfinite(reported2[base]):
        assert np.isclose(reported[base], reported2[base]), (
            f"{base} changed when QC annotations were removed: "
            f"{reported[base]} vs {reported2[base]}"
        )

print(">> All checks passed", flush=True)
