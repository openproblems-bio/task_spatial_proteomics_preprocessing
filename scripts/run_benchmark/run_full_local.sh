#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# NOTE: depending on the the datasets and components, you may need to launch this workflow
# on a different compute platform (e.g. a HPC, AWS Cloud, Azure Cloud, Google Cloud).
# please refer to the nextflow information for more details:
# https://www.nextflow.io/docs/latest/

set -e

echo "Running benchmark"
echo "  Make sure to run 'scripts/project/build_all_docker_containers.sh'!"

# The dataset to benchmark against: a SpatialData zarr with a multiplexed image
# and, optionally, reference labels. `scripts/create_resources/isp/prepare_exemplar_001.sh`
# builds one from the public MCMICRO exemplar-001 acquisition.
INPUT_ZARR="${INPUT_ZARR:?Set INPUT_ZARR=/path/to/dataset.zarr}"
RUN_NAME="${RUN_NAME:-exemplar_001}"

# generate a unique id
RUN_ID="run_$(date +%Y-%m-%d_%H-%M-%S)"
publish_dir="resources/results/${RUN_ID}"

# write the parameters to file
cat > /tmp/params.yaml << HERE
id: "$RUN_NAME"
input_zarr: "$INPUT_ZARR"
output_state: "state.yaml"
publish_dir: "$publish_dir"
HERE

# Panel-dependent metric arguments have no safe default, so append them here for
# the panel being run. For MCMICRO exemplar-001:
#
#   nuclear_markers: "DNA_6,DNA_7,DNA_8"
#   exclusive_pairs: "ECAD:CD45,ECAD:CD11B,SMA:CD45"
#
# Without them the nuclear-agreement and co-expression metrics report NaN rather
# than guessing for an unknown panel.
if [ -n "${NUCLEAR_MARKERS:-}" ]; then
  echo "nuclear_markers: \"$NUCLEAR_MARKERS\"" >> /tmp/params.yaml
fi
if [ -n "${EXCLUSIVE_PAIRS:-}" ]; then
  echo "exclusive_pairs: \"$EXCLUSIVE_PAIRS\"" >> /tmp/params.yaml
fi
if [ -n "${METHOD_PARAMETERS_YAML:-}" ]; then
  echo "method_parameters_yaml: \"$METHOD_PARAMETERS_YAML\"" >> /tmp/params.yaml
fi

# run the benchmark
nextflow run openproblems-bio/task_spatial_proteomics_preprocessing \
  -r build/main \
  -main-script target/nextflow/workflows/run_benchmark/main.nf \
  -profile docker \
  -resume \
  -c common/nextflow_helpers/labels_ci.config \
  -params-file /tmp/params.yaml
