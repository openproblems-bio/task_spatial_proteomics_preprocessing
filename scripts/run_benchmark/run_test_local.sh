#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

set -e

# A small SpatialData zarr to smoke-test the workflow against. There is no
# `resources_test` bucket for this task yet, so build one locally first:
#
#   bash scripts/create_resources/isp/prepare_exemplar_001.sh \
#     resources_test/task_spatial_proteomics_preprocessing/exemplar_001
#
INPUT_ZARR="${INPUT_ZARR:-resources_test/task_spatial_proteomics_preprocessing/exemplar_001/dataset.zarr}"

if [ ! -d "$INPUT_ZARR" ]; then
  echo "ERROR: no input zarr at $INPUT_ZARR" >&2
  echo "Build one with scripts/create_resources/isp/prepare_exemplar_001.sh," >&2
  echo "or set INPUT_ZARR to an existing SpatialData store." >&2
  exit 1
fi

echo "Running benchmark on test data"
echo "  Make sure to run 'scripts/project/build_all_docker_containers.sh'!"

# generate a unique id
RUN_ID="testrun_$(date +%Y-%m-%d_%H-%M-%S)"
publish_dir="temp/results/${RUN_ID}"

# Only the baseline chain and one segmentation alternative, to keep the smoke
# test short; the full run is scripts/run_benchmark/run_full_local.sh.
nextflow run . \
  -main-script target/nextflow/workflows/run_benchmark/main.nf \
  -profile docker \
  -resume \
  -c common/nextflow_helpers/labels_ci.config \
  --id exemplar_001 \
  --input_zarr "$INPUT_ZARR" \
  --image_correction_methods basic_image_correction \
  --segmentation_methods "basic_threshold_segmentation;watershed_expand_segmentation" \
  --aggregation_methods mean_intensity_aggregation \
  --qc_methods basic_qc_filter \
  --normalization_methods total_log1p_normalization \
  --output_state state.yaml \
  --publish_dir "$publish_dir"
