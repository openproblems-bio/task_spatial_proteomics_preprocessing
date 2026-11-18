#!/usr/bin/env bash

set -euo pipefail

REPO_ROOT="$(git rev-parse --show-toplevel 2>/dev/null || pwd)"
cd "$REPO_ROOT"

OUTPUT_DIR="${OUTPUT_DIR:-output/temp_isp}"
ZARR_PATH="${ZARR_PATH:?Please set ZARR_PATH=/path/to/exemplar-001.zarr}"
IMAGE_KEY="${IMAGE_KEY:-null}"
IMAGE_CORRECTION_METHOD="${IMAGE_CORRECTION_METHOD:-percentile_background}"
SEGMENTATION_METHOD="${SEGMENTATION_METHOD:-threshold}"
AGGREGATION_METHOD="${AGGREGATION_METHOD:-mean}"
QC_METHOD="${QC_METHOD:-fixed_threshold}"
NORMALIZATION_METHOD="${NORMALIZATION_METHOD:-total_log1p}"
LABELS_KEY="${LABELS_KEY:-null}"
# Components must be built with the same engine the workflow runs them under.
VIASH_ENGINE="${VIASH_ENGINE:-docker}"

mkdir -p "$OUTPUT_DIR"

echo "Using repository:              $REPO_ROOT"
echo "Using Zarr input:              $ZARR_PATH"
echo "Using output dir:              $OUTPUT_DIR"
echo "Using image correction method: $IMAGE_CORRECTION_METHOD"
echo "Using segmentation method:     $SEGMENTATION_METHOD"
echo "Using aggregation method:      $AGGREGATION_METHOD"
echo "Using QC method:               $QC_METHOD"
echo "Using normalization method:    $NORMALIZATION_METHOD"
echo "Using viash engine:            $VIASH_ENGINE"

if [ ! -d "$ZARR_PATH" ]; then
  echo "ERROR: ZARR_PATH does not exist or is not a directory: $ZARR_PATH" >&2
  echo "Set it explicitly, for example:" >&2
  echo "ZARR_PATH=/path/to/exemplar-001.zarr bash src/workflows/run_benchmark/test.sh" >&2
  exit 1
fi

echo "Building required Viash components..."

for q in \
  "methods_image_correction/no_image_correction" \
  "methods_image_correction/basic_image_correction" \
  "methods_image_correction/flatfield_image_correction" \
  "methods_proteomics_segmentation/basic_threshold_segmentation" \
  "methods_proteomics_segmentation/watershed_expand_segmentation" \
  "methods_proteomics_segmentation/reference_segmentation" \
  "methods_protein_aggregation/mean_intensity_aggregation" \
  "methods_protein_aggregation/sum_intensity_aggregation" \
  "methods_protein_aggregation/median_intensity_aggregation" \
  "methods_proteomics_qc/no_qc_filter" \
  "methods_proteomics_qc/basic_qc_filter" \
  "methods_proteomics_qc/adaptive_quantile_qc_filter" \
  "methods_proteomics_normalization/no_normalization" \
  "methods_proteomics_normalization/total_log1p_normalization" \
  "methods_proteomics_normalization/arcsinh_normalization" \
  "metrics_proteomics/basic" \
  "metrics_proteomics/marker_distribution" \
  "metrics_proteomics/normalization_diagnostics"  \
  "metrics_proteomics/segmentation_quality" \
  "metrics_proteomics/technical_covariates" \
  "metrics_proteomics/signal_background"
do
  viash ns build --engine "$VIASH_ENGINE" --query "$q"
done

echo "Running ISP preprocessing benchmark workflow..."

export NXF_VER="${NXF_VER:-24.04.3}"

nextflow run src/workflows/run_benchmark/main.nf \
  -entry run_wf \
  -resume \
  --zarr "$ZARR_PATH" \
  --outdir "$OUTPUT_DIR" \
  --viash_engine "$VIASH_ENGINE" \
  --image_key "$IMAGE_KEY" \
  --image_correction_method "$IMAGE_CORRECTION_METHOD" \
  --background_percentile "${BACKGROUND_PERCENTILE:-1.0}" \
  --rescale_max "${RESCALE_MAX:-false}" \
  --flatfield_sigma "${FLATFIELD_SIGMA:-50.0}" \
  --flatfield_mode "${FLATFIELD_MODE:-subtract}" \
  --segmentation_method "$SEGMENTATION_METHOD" \
  --labels_key "$LABELS_KEY" \
  --nuclei_channel "${NUCLEI_CHANNEL:-0}" \
  --threshold_method "${THRESHOLD_METHOD:-quantile}" \
  --threshold_quantile "${THRESHOLD_QUANTILE:-0.90}" \
  --smoothing_sigma "${SMOOTHING_SIGMA:-1.0}" \
  --segmentation_min_area "${SEGMENTATION_MIN_AREA:-30}" \
  --hole_area "${HOLE_AREA:-64}" \
  --expansion_distance "${EXPANSION_DISTANCE:-5}" \
  --aggregation_method "$AGGREGATION_METHOD" \
  --stat "${STAT:-mean}" \
  --qc_method "$QC_METHOD" \
  --min_area "${MIN_AREA:-50}" \
  --min_total "${MIN_TOTAL:-0.0}" \
  --area_quantile "${AREA_QUANTILE:-0.05}" \
  --total_quantile "${TOTAL_QUANTILE:-0.05}" \
  --normalization_method "$NORMALIZATION_METHOD" \
  --scale "${SCALE:-10000}" \
  --cofactor "${COFACTOR:-5.0}" \
  --transform "${TRANSFORM:-log1p}" \
  --cov_thr "${COV_THR:-0.1}" \
  --metric_marker_layer "${METRIC_MARKER_LAYER:-counts}" \
  --zero_threshold "${ZERO_THRESHOLD:-0.0}" \
  --metric_low_quantile "${METRIC_LOW_QUANTILE:-0.01}" \
  --metric_high_quantile "${METRIC_HIGH_QUANTILE:-0.99}" \
  --metric_raw_layer "${METRIC_RAW_LAYER:-counts}" \
  --metric_normalized_layer "${METRIC_NORMALIZED_LAYER:-normalized}" \
  --metric_max_cells "${METRIC_MAX_CELLS:-5000}" \
  --nuclear_markers "${NUCLEAR_MARKERS:-null}" \
  --exclusive_pairs "${EXCLUSIVE_PAIRS:-null}" \
  --n_neighbors "${N_NEIGHBORS:-6}" \
  --metric_covariate "${METRIC_COVARIATE:-area}" \
  --metric_n_pcs "${METRIC_N_PCS:-10}" \
  --dilation_iterations "${DILATION_ITERATIONS:-3}" \
  --norm_percentile "${NORM_PERCENTILE:-95.0}" \
  --metric_random_state "${METRIC_RANDOM_STATE:-0}"

echo
echo "Workflow finished. Published outputs:"
ls -lh "$OUTPUT_DIR"
