#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# This needs s3://openproblems-data/resources_test/task_spatial_proteomics_preprocessing/,
# which has not been populated yet - nothing has been uploaded for this task. Until
# it is, run scripts/run_benchmark/run_test_local.sh against a locally built zarr.
echo "TODO: upload test resources before using this script."
echo "  See the test_resources note in _viash.yaml."
exit 1

set -e

resources_test_s3=s3://openproblems-data/resources_test/task_spatial_proteomics_preprocessing
publish_dir_s3="s3://openproblems-nextflow/temp/results/$(date +%Y-%m-%d_%H-%M-%S)"

# write the parameters to file
cat > /tmp/params.yaml << HERE
id: exemplar_001
input_zarr: $resources_test_s3/exemplar_001/dataset.zarr
output_state: "state.yaml"
publish_dir: $publish_dir_s3
HERE

tw launch https://github.com/openproblems-bio/task_spatial_proteomics_preprocessing.git \
  --revision build/main \
  --pull-latest \
  --main-script target/nextflow/workflows/run_benchmark/main.nf \
  --workspace 53907369739130 \
  --compute-env 6TeIFgV5OY4pJCk8I0bfOh \
  --params-file /tmp/params.yaml \
  --config common/nextflow_helpers/labels_tw.config \
  --labels task_spatial_proteomics_preprocessing,test
