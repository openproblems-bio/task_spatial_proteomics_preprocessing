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

# Methods with the classical threshold arguments; only these accept them.
THRESHOLD_METHODS = {"basic_threshold_segmentation", "watershed_expand_segmentation"}
# Per-method arguments that fit the 64x64 fixture, whose nuclei are ~10 px across.
FIXTURE_ARGS = {
    # min_area is set above the smallest fixture blob so that object removal,
    # and therefore the label-renumbering path, is actually exercised.
    "basic_threshold_segmentation": ["--min_area", "20"],
    "watershed_expand_segmentation": ["--min_area", "20"],
    "watershed_txsim_segmentation": ["--min_size", "20", "--threshold_func", "otsu"],
    "binning_segmentation": ["--bin_size", "16"],
    "cellpose_segmentation": ["--diameter", "10", "--model_type", "nuclei"],
    "cellposev4_segmentation": ["--diameter", "10"],
    "stardist_segmentation": [],
}

if name == "reference_segmentation":
    zarr_path = work / "fixture.zarr"
    isp_fixtures.make_spatialdata_zarr(zarr_path)
    cmd += ["--input_zarr", str(zarr_path), "--labels_key", "fixture_cell"]
else:
    image_path = work / "corrected_image.npy"
    markers_path = work / "corrected_markers.json"
    isp_fixtures.make_image(image_path)
    isp_fixtures.make_markers_json(markers_path)
    cmd += ["--input_image", str(image_path), "--nuclei_channel", "0"] + FIXTURE_ARGS.get(name, [])

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

if name == "binning_segmentation":
    print(">> Checking the bin grid", flush=True)
    n_bins = -(-labels.shape[0] // 16) * -(-labels.shape[1] // 16)
    assert present.size == n_bins, f"expected {n_bins} bins, got {present.size}"


def run_variant(tag, extra, overrides=None):
    """Re-run the component with extra arguments, writing to separate outputs.

    `overrides` replaces the value of flags already in `cmd` (e.g. a different
    --model_type), rather than appending a second copy of the flag.
    """
    out_labels = work / f"labels_{tag}.npy"
    out_meta = work / f"metadata_{tag}.json"
    base = list(cmd)
    base[base.index("--output_labels") + 1] = str(out_labels)
    base[base.index("--output_metadata") + 1] = str(out_meta)
    for flag, value in (overrides or {}).items():
        if flag in base:
            base[base.index(flag) + 1] = value
        else:
            base += [flag, value]
    res = subprocess.run(base + extra, capture_output=True, text=True)
    if res.returncode != 0:
        return res, None, None
    with out_meta.open() as f:
        return res, np.load(out_labels), json.load(f)


if name != "reference_segmentation" and name != "binning_segmentation":
    print(">> Selecting the nuclear channel by marker name", flush=True)
    res, by_name, meta_name = run_variant("by_name", ["--input_markers", str(markers_path), "--nuclear_markers", "DNA_1"])
    assert res.returncode == 0, res.stderr
    # DNA_1 is channel 0, so naming it must give the same mask as the index.
    assert np.array_equal(by_name, labels), "--nuclear_markers DNA_1 differs from --nuclei_channel 0"
    assert meta_name["input_channels"]["nuclear"]["markers"] == ["DNA_1"]

    print(">> Max-projecting several nuclear markers", flush=True)
    res, two, meta_two = run_variant("two", ["--input_markers", str(markers_path), "--nuclear_markers", "dna_1, DNA_2"])
    assert res.returncode == 0, res.stderr
    assert two.shape == labels.shape
    # Both markers must actually be used, resolved case-insensitively.
    assert meta_two["input_channels"]["nuclear"] == {"markers": ["DNA_1", "DNA_2"], "channels": [0, 1]}, (
        meta_two["input_channels"]
    )
    assert meta_two["n_objects"] > 0, "max-projecting two replicate nuclear stains lost every object"

    print(">> Rejecting an unknown marker name", flush=True)
    res, _, _ = run_variant("bad", ["--input_markers", str(markers_path), "--nuclear_markers", "NOT_A_MARKER"])
    assert res.returncode != 0, "an unknown marker name must fail instead of falling back to the index"

    print(">> Rejecting marker names without a marker sidecar", flush=True)
    res, _, _ = run_variant("no_sidecar", ["--nuclear_markers", "DNA_1"])
    assert res.returncode != 0, "marker names without --input_markers must fail"

if name in ("cellpose_segmentation", "cellposev4_segmentation"):
    print(">> Running with a membrane channel", flush=True)
    membrane_args = ["--input_markers", str(markers_path), "--nuclear_markers", "DNA_1", "--membrane_markers", "ECAD"]
    # The fixture's nuclei model cannot use a membrane channel, so use a cyto model.
    overrides = {"--model_type": "cyto"} if name == "cellpose_segmentation" else None
    res, mem_labels, meta_mem = run_variant("membrane", membrane_args, overrides)
    assert res.returncode == 0, res.stderr
    assert mem_labels.shape == labels.shape
    assert meta_mem["input_channels"]["membrane"] == {"markers": ["ECAD"], "channels": [3]}, meta_mem["input_channels"]
    if name == "cellpose_segmentation":
        assert meta_mem["channels"] == [1, 2], "two-channel mode must be [membrane, nuclear]"

        print(">> The nuclei model ignores membrane markers", flush=True)
        res, _, meta_nuc = run_variant("nuclei_membrane", membrane_args)
        assert res.returncode == 0, res.stderr
        assert meta_nuc["input_channels"]["membrane"] is None
        assert meta_nuc["channels"] == [0, 0]

if name in THRESHOLD_METHODS:
    print(">> Rejecting an invalid --threshold_method", flush=True)
    bad = subprocess.run(
        cmd + ["--threshold_method", "not_a_method"],
        capture_output=True,
        text=True,
    )
    assert bad.returncode != 0, "an unknown --threshold_method must fail instead of silently falling back"

print(">> All checks passed", flush=True)
