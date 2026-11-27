#!/usr/bin/env bash
#
# Run the benchmark on a local SpatialData zarr.
#
#   ZARR_PATH=/path/to/exemplar-001.zarr \
#   NUCLEAR_MARKERS="DNA_6,DNA_7,DNA_8" \
#   EXCLUSIVE_PAIRS="ECAD:CD45,ECAD:CD11B,SMA:CD45" \
#   bash src/workflows/run_benchmark/test.sh
#
# The method combinations are chosen by the workflow, not here: with the shipped
# defaults this runs the baseline chain plus one variant per alternative method.
# To run a different set, pass e.g.
#
#   SEGMENTATION_METHODS="basic_threshold_segmentation;watershed_expand_segmentation"
#
# and to vary a method's parameters, point METHOD_PARAMETERS_YAML at a file.

set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
cd "$REPO_ROOT"

ZARR_PATH="${ZARR_PATH:?Please set ZARR_PATH=/path/to/exemplar-001.zarr}"
OUTPUT_DIR="${OUTPUT_DIR:-output/temp_isp}"
RUN_ID="${RUN_ID:-exemplar_001}"

# Nextflow profile selecting how components are executed. The workflow no longer
# shells out to the built executables, so the engine is chosen here rather than
# with a workflow argument: use `singularity` where Docker is unavailable.
PROFILE="${PROFILE:-docker}"

if [ ! -d "$ZARR_PATH" ]; then
  echo "ERROR: ZARR_PATH does not exist or is not a directory: $ZARR_PATH" >&2
  exit 1
fi

echo "Using repository:   $REPO_ROOT"
echo "Using Zarr input:   $ZARR_PATH"
echo "Using output dir:   $OUTPUT_DIR"
echo "Using profile:      $PROFILE"

echo "Building components and workflow..."
if [ "$PROFILE" = "docker" ]; then
  # also build the container images the docker profile will run
  viash ns build --setup cachedbuild --parallel
else
  viash ns build --parallel
fi

args=(
  -main-script target/nextflow/workflows/run_benchmark/main.nf
  -profile "$PROFILE"
  -resume
  -c common/nextflow_helpers/labels_ci.config
  --id "$RUN_ID"
  --input_zarr "$ZARR_PATH"
  --output_state state.yaml
  --publish_dir "$OUTPUT_DIR"
)

# Optional arguments, passed only when set, so the workflow's own defaults apply.
# `add_opt NAME --flag` appends nothing when $NAME is unset or empty.
add_opt() {
  local value="${!1:-}"
  if [ -n "$value" ]; then
    args+=("$2" "$value")
  fi
}

add_opt IMAGE_KEY                --image_key
add_opt LABELS_KEY               --labels_key
add_opt NUCLEAR_MARKERS          --nuclear_markers
add_opt EXCLUSIVE_PAIRS          --exclusive_pairs
add_opt DILATION_ITERATIONS      --dilation_iterations
add_opt METHOD_PARAMETERS_YAML   --method_parameters_yaml
add_opt DEFAULT_METHODS          --default_methods
add_opt IMAGE_CORRECTION_METHODS --image_correction_methods
add_opt SEGMENTATION_METHODS     --segmentation_methods
add_opt AGGREGATION_METHODS      --aggregation_methods
add_opt QC_METHODS               --qc_methods
add_opt NORMALIZATION_METHODS    --normalization_methods

echo "Running spatial proteomics preprocessing benchmark..."
export NXF_VER="${NXF_VER:-24.04.3}"
nextflow run . "${args[@]}"

echo
echo "Workflow finished. Published outputs:"
ls -lh "$OUTPUT_DIR"
