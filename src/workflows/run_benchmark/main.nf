nextflow.enable.dsl = 2

/*
 * Spatial proteomics preprocessing benchmark workflow.
 *
 * Current modular workflow:
 * SpatialData zarr
 *   -> image correction
 *   -> segmentation
 *   -> cell-by-protein aggregation
 *   -> quality control
 *   -> normalization
 *   -> preprocessing metrics
 *
 * The workflow supports benchmarkable method choices for image correction,
 * segmentation, protein aggregation, quality control, and normalization using
 * --image_correction_method, --segmentation_method, --aggregation_method,
 * --qc_method, and --normalization_method. Metrics are computed as three
 * complementary benchmark families: basic summary metrics, marker distribution
 * metrics, and normalization diagnostics.
 */

ROOT = workflow.launchDir
PUBLISH_DIR = "${params.outdir ?: params.publish_dir ?: 'results_isp'}"

// Engine the Viash executables run under. Docker is the default and what the
// benchmark is normally run with; 'native' lets the workflow run where Docker
// is unavailable, provided the components were built with the same engine.
VIASH_ENGINE = "${params.viash_engine ?: 'docker'}"

NO_IMAGE_CORRECTION_EXE = "${ROOT}/target/executable/methods_image_correction/no_image_correction/no_image_correction"
BASIC_IMAGE_CORRECTION_EXE = "${ROOT}/target/executable/methods_image_correction/basic_image_correction/basic_image_correction"
FLATFIELD_IMAGE_CORRECTION_EXE = "${ROOT}/target/executable/methods_image_correction/flatfield_image_correction/flatfield_image_correction"

THRESHOLD_SEGMENTATION_EXE = "${ROOT}/target/executable/methods_proteomics_segmentation/basic_threshold_segmentation/basic_threshold_segmentation"
WATERSHED_SEGMENTATION_EXE = "${ROOT}/target/executable/methods_proteomics_segmentation/watershed_expand_segmentation/watershed_expand_segmentation"
REFERENCE_SEGMENTATION_EXE = "${ROOT}/target/executable/methods_proteomics_segmentation/reference_segmentation/reference_segmentation"

MEAN_AGGREGATION_EXE = "${ROOT}/target/executable/methods_protein_aggregation/mean_intensity_aggregation/mean_intensity_aggregation"
SUM_AGGREGATION_EXE = "${ROOT}/target/executable/methods_protein_aggregation/sum_intensity_aggregation/sum_intensity_aggregation"
MEDIAN_AGGREGATION_EXE = "${ROOT}/target/executable/methods_protein_aggregation/median_intensity_aggregation/median_intensity_aggregation"

NO_QC_EXE = "${ROOT}/target/executable/methods_proteomics_qc/no_qc_filter/no_qc_filter"
BASIC_QC_EXE = "${ROOT}/target/executable/methods_proteomics_qc/basic_qc_filter/basic_qc_filter"
ADAPTIVE_QC_EXE = "${ROOT}/target/executable/methods_proteomics_qc/adaptive_quantile_qc_filter/adaptive_quantile_qc_filter"
NO_NORMALIZATION_EXE = "${ROOT}/target/executable/methods_proteomics_normalization/no_normalization/no_normalization"
TOTAL_LOG1P_NORMALIZATION_EXE = "${ROOT}/target/executable/methods_proteomics_normalization/total_log1p_normalization/total_log1p_normalization"
ARCSINH_NORMALIZATION_EXE = "${ROOT}/target/executable/methods_proteomics_normalization/arcsinh_normalization/arcsinh_normalization"
METRICS_BASIC_EXE = "${ROOT}/target/executable/metrics_proteomics/basic/basic"
METRICS_MARKER_DISTRIBUTION_EXE = "${ROOT}/target/executable/metrics_proteomics/marker_distribution/marker_distribution"
METRICS_NORMALIZATION_DIAGNOSTICS_EXE = "${ROOT}/target/executable/metrics_proteomics/normalization_diagnostics/normalization_diagnostics"
METRICS_SEGMENTATION_QUALITY_EXE = "${ROOT}/target/executable/metrics_proteomics/segmentation_quality/segmentation_quality"
METRICS_TECHNICAL_COVARIATES_EXE = "${ROOT}/target/executable/metrics_proteomics/technical_covariates/technical_covariates"
METRICS_SIGNAL_BACKGROUND_EXE = "${ROOT}/target/executable/metrics_proteomics/signal_background/signal_background"

process image_correction {
    tag "image_correction_${params.image_correction_method ?: 'percentile_background'}"
    publishDir PUBLISH_DIR, mode: 'copy'

    input:
    path zarr

    output:
    path "corrected_image.npy", emit: image
    path "corrected_markers.json", emit: markers

    script:
    def method = params.image_correction_method ?: 'percentile_background'
    def image_key = params.image_key ?: 'null'
    def background_percentile = params.background_percentile ?: 1.0
    def rescale_max = params.rescale_max ?: false
    def flatfield_sigma = params.flatfield_sigma ?: 50.0
    def flatfield_mode = params.flatfield_mode ?: 'subtract'
    """
    set -euo pipefail
    export VIASH_ENGINE_TYPE=${VIASH_ENGINE}

    METHOD="${method}"
    ZARR_PATH="\$(realpath '${zarr}')"
    OUT_IMAGE="\$PWD/corrected_image.npy"
    OUT_MARKERS="\$PWD/corrected_markers.json"

    case "\$METHOD" in
      no_correction|none)
        "${NO_IMAGE_CORRECTION_EXE}" \
          --input_zarr "\$ZARR_PATH" \
          --image_key "${image_key}" \
          --output_image "\$OUT_IMAGE" \
          --output_markers "\$OUT_MARKERS"
        ;;
      percentile_background|basic)
        "${BASIC_IMAGE_CORRECTION_EXE}" \
          --input_zarr "\$ZARR_PATH" \
          --image_key "${image_key}" \
          --background_percentile "${background_percentile}" \
          --rescale_max "${rescale_max}" \
          --output_image "\$OUT_IMAGE" \
          --output_markers "\$OUT_MARKERS"
        ;;
      flatfield|flatfield_correction)
        "${FLATFIELD_IMAGE_CORRECTION_EXE}" \
          --input_zarr "\$ZARR_PATH" \
          --image_key "${image_key}" \
          --background_percentile "${background_percentile}" \
          --rescale_max "${rescale_max}" \
          --flatfield_sigma "${flatfield_sigma}" \
          --flatfield_mode "${flatfield_mode}" \
          --output_image "\$OUT_IMAGE" \
          --output_markers "\$OUT_MARKERS"
        ;;
      *)
        echo "ERROR: Unknown image correction method: \$METHOD" >&2
        echo "Use one of: no_correction, percentile_background, flatfield" >&2
        exit 1
        ;;
    esac
    """
}

process segment_cells {
    tag "segmentation_${params.segmentation_method ?: 'threshold'}"
    publishDir PUBLISH_DIR, mode: 'copy'

    input:
    path corrected_image
    path zarr

    output:
    path "segmentation_labels.npy", emit: labels
    path "segmentation_metadata.json", optional: true, emit: metadata

    script:
    def method = params.segmentation_method ?: 'threshold'
    def nuclei_channel = params.nuclei_channel ?: 0
    def threshold_quantile = params.threshold_quantile ?: 0.90
    def min_area = params.segmentation_min_area ?: 30
    def threshold_method = params.threshold_method ?: 'quantile'
    def smoothing_sigma = params.smoothing_sigma ?: 1.0
    def hole_area = params.hole_area ?: 64
    def expansion_distance = params.expansion_distance ?: 5
    def labels_key = params.labels_key ?: 'null'
    """
    set -euo pipefail
    export VIASH_ENGINE_TYPE=${VIASH_ENGINE}

    METHOD="${method}"
    IMAGE_PATH="\$(realpath '${corrected_image}')"
    ZARR_PATH="\$(realpath '${zarr}')"
    OUT_LABELS="\$PWD/segmentation_labels.npy"

    case "\$METHOD" in
      threshold|basic_threshold)
        "${THRESHOLD_SEGMENTATION_EXE}" \
          --input_image "\$IMAGE_PATH" \
          --nuclei_channel "${nuclei_channel}" \
          --threshold_method "${threshold_method}" \
          --threshold_quantile "${threshold_quantile}" \
          --min_area "${min_area}" \
          --output_labels "\$OUT_LABELS" \
          --output_metadata "\$PWD/segmentation_metadata.json"
        ;;
      watershed_expand|watershed)
        "${WATERSHED_SEGMENTATION_EXE}" \
          --input_image "\$IMAGE_PATH" \
          --nuclei_channel "${nuclei_channel}" \
          --threshold_method "${threshold_method}" \
          --threshold_quantile "${threshold_quantile}" \
          --smoothing_sigma "${smoothing_sigma}" \
          --min_area "${min_area}" \
          --hole_area "${hole_area}" \
          --expansion_distance "${expansion_distance}" \
          --output_labels "\$OUT_LABELS" \
          --output_metadata "\$PWD/segmentation_metadata.json"
        ;;
      reference|existing_labels)
        "${REFERENCE_SEGMENTATION_EXE}" \
          --input_zarr "\$ZARR_PATH" \
          --labels_key "${labels_key}" \
          --output_labels "\$OUT_LABELS" \
          --output_metadata "\$PWD/segmentation_metadata.json"
        ;;
      *)
        echo "ERROR: Unknown segmentation method: \$METHOD" >&2
        echo "Use one of: threshold, watershed_expand, reference" >&2
        exit 1
        ;;
    esac
    """
}

process aggregate_proteins {
    tag "aggregation_${params.aggregation_method ?: 'mean'}"
    publishDir PUBLISH_DIR, mode: 'copy'

    input:
    path corrected_image
    path segmentation_labels
    path corrected_markers

    output:
    path "cell_by_protein.h5ad", emit: adata

    script:
    def method = params.aggregation_method ?: 'mean'
    """
    set -euo pipefail
    export VIASH_ENGINE_TYPE=${VIASH_ENGINE}

    METHOD="${method}"
    IMAGE_PATH="\$(realpath '${corrected_image}')"
    LABELS_PATH="\$(realpath '${segmentation_labels}')"
    MARKERS_PATH="\$(realpath '${corrected_markers}')"
    OUT_PATH="\$PWD/cell_by_protein.h5ad"

    case "\$METHOD" in
      mean|mean_intensity)
        "${MEAN_AGGREGATION_EXE}" \
          --input_image "\$IMAGE_PATH" \
          --input_labels "\$LABELS_PATH" \
          --input_markers "\$MARKERS_PATH" \
          --output "\$OUT_PATH"
        ;;
      sum|sum_intensity)
        "${SUM_AGGREGATION_EXE}" \
          --input_image "\$IMAGE_PATH" \
          --input_labels "\$LABELS_PATH" \
          --input_markers "\$MARKERS_PATH" \
          --output "\$OUT_PATH"
        ;;
      median|median_intensity)
        "${MEDIAN_AGGREGATION_EXE}" \
          --input_image "\$IMAGE_PATH" \
          --input_labels "\$LABELS_PATH" \
          --input_markers "\$MARKERS_PATH" \
          --output "\$OUT_PATH"
        ;;
      *)
        echo "ERROR: Unknown aggregation method: \$METHOD" >&2
        echo "Use one of: mean, sum, median" >&2
        exit 1
        ;;
    esac
    """
}

process qc_filter {
    tag "qc_${params.qc_method ?: 'fixed_threshold'}"
    publishDir PUBLISH_DIR, mode: 'copy'

    input:
    path adata_in

    output:
    path "cell_by_protein_qc.h5ad", emit: adata

    script:
    def method = params.qc_method ?: 'fixed_threshold'
    def min_area = params.containsKey('min_area') ? params.min_area : 50
    def min_total = params.containsKey('min_total') ? params.min_total : 0.0
    def area_quantile = params.containsKey('area_quantile') ? params.area_quantile : 0.05
    def total_quantile = params.containsKey('total_quantile') ? params.total_quantile : 0.05
    """
    set -euo pipefail
    export VIASH_ENGINE_TYPE=${VIASH_ENGINE}

    METHOD="${method}"
    IN_PATH="\$(realpath '${adata_in}')"
    OUT_PATH="\$PWD/cell_by_protein_qc.h5ad"

    case "\$METHOD" in
      no_qc|none|all_pass)
        "${NO_QC_EXE}" \
          --input "\$IN_PATH" \
          --min_area "${min_area}" \
          --min_total "${min_total}" \
          --area_quantile "${area_quantile}" \
          --total_quantile "${total_quantile}" \
          --output "\$OUT_PATH"
        ;;
      fixed_threshold|basic)
        "${BASIC_QC_EXE}" \
          --input "\$IN_PATH" \
          --min_area "${min_area}" \
          --min_total "${min_total}" \
          --area_quantile "${area_quantile}" \
          --total_quantile "${total_quantile}" \
          --output "\$OUT_PATH"
        ;;
      adaptive_quantile|quantile)
        "${ADAPTIVE_QC_EXE}" \
          --input "\$IN_PATH" \
          --min_area "${min_area}" \
          --min_total "${min_total}" \
          --area_quantile "${area_quantile}" \
          --total_quantile "${total_quantile}" \
          --output "\$OUT_PATH"
        ;;
      *)
        echo "ERROR: Unknown QC method: \$METHOD" >&2
        echo "Use one of: no_qc, fixed_threshold, adaptive_quantile" >&2
        exit 1
        ;;
    esac
    """
}

process normalize {
    tag "normalization_${params.normalization_method ?: 'total_log1p'}"
    publishDir PUBLISH_DIR, mode: 'copy'

    input:
    path adata_in

    output:
    path "cell_by_protein_norm.h5ad", emit: adata

    script:
    def method = params.normalization_method ?: 'total_log1p'
    def scale = params.containsKey('scale') ? params.scale : 10000
    def cofactor = params.containsKey('cofactor') ? params.cofactor : 5.0
    """
    set -euo pipefail
    export VIASH_ENGINE_TYPE=${VIASH_ENGINE}

    METHOD="${method}"
    IN_PATH="\$(realpath '${adata_in}')"
    OUT_PATH="\$PWD/cell_by_protein_norm.h5ad"

    case "\$METHOD" in
      no_normalization|none|raw)
        "${NO_NORMALIZATION_EXE}" \
          --input "\$IN_PATH" \
          --scale "${scale}" \
          --cofactor "${cofactor}" \
          --output "\$OUT_PATH"
        ;;
      total_log1p|total_log1p_normalization|log1p|total)
        "${TOTAL_LOG1P_NORMALIZATION_EXE}" \
          --input "\$IN_PATH" \
          --scale "${scale}" \
          --cofactor "${cofactor}" \
          --output "\$OUT_PATH"
        ;;
      arcsinh|arcsinh_normalization|asinh)
        "${ARCSINH_NORMALIZATION_EXE}" \
          --input "\$IN_PATH" \
          --scale "${scale}" \
          --cofactor "${cofactor}" \
          --output "\$OUT_PATH"
        ;;
      *)
        echo "ERROR: Unknown normalization method: \$METHOD" >&2
        echo "Use one of: no_normalization, total_log1p, arcsinh" >&2
        exit 1
        ;;
    esac
    """
}

process compute_basic_metrics {
    tag "metrics_basic"
    publishDir PUBLISH_DIR, mode: 'copy'

    input:
    path adata_in

    output:
    path "metrics_basic.h5ad", emit: metrics

    script:
    def cov_thr = params.cov_thr ?: 0.1
    """
    set -euo pipefail
    export VIASH_ENGINE_TYPE=${VIASH_ENGINE}

    IN_PATH="\$(realpath '${adata_in}')"
    OUT_PATH="\$PWD/metrics_basic.h5ad"

    "${METRICS_BASIC_EXE}" \
      --input "\$IN_PATH" \
      --cov_thr "${cov_thr}" \
      --output "\$OUT_PATH"
    """
}

process compute_marker_distribution_metrics {
    tag "metrics_marker_distribution"
    publishDir PUBLISH_DIR, mode: 'copy'

    input:
    path adata_in

    output:
    path "metrics_marker_distribution.h5ad", emit: metrics

    script:
    def layer = params.metric_marker_layer ?: 'counts'
    def zero_threshold = params.containsKey('zero_threshold') ? params.zero_threshold : 0.0
    def coverage_threshold = params.cov_thr ?: 0.1
    def low_quantile = params.containsKey('metric_low_quantile') ? params.metric_low_quantile : 0.01
    def high_quantile = params.containsKey('metric_high_quantile') ? params.metric_high_quantile : 0.99
    """
    set -euo pipefail
    export VIASH_ENGINE_TYPE=${VIASH_ENGINE}

    IN_PATH="\$(realpath '${adata_in}')"
    OUT_PATH="\$PWD/metrics_marker_distribution.h5ad"

    "${METRICS_MARKER_DISTRIBUTION_EXE}" \
      --input "\$IN_PATH" \
      --layer "${layer}" \
      --zero_threshold "${zero_threshold}" \
      --coverage_threshold "${coverage_threshold}" \
      --low_quantile "${low_quantile}" \
      --high_quantile "${high_quantile}" \
      --output "\$OUT_PATH"
    """
}

process compute_normalization_diagnostics {
    tag "metrics_normalization_diagnostics"
    publishDir PUBLISH_DIR, mode: 'copy'

    input:
    path adata_in

    output:
    path "metrics_normalization_diagnostics.h5ad", emit: metrics

    script:
    def raw_layer = params.metric_raw_layer ?: 'counts'
    def normalized_layer = params.metric_normalized_layer ?: 'normalized'
    def max_cells = params.containsKey('metric_max_cells') ? params.metric_max_cells : 5000
    def random_state = params.containsKey('metric_random_state') ? params.metric_random_state : 0
    """
    set -euo pipefail
    export VIASH_ENGINE_TYPE=${VIASH_ENGINE}

    IN_PATH="\$(realpath '${adata_in}')"
    OUT_PATH="\$PWD/metrics_normalization_diagnostics.h5ad"

    "${METRICS_NORMALIZATION_DIAGNOSTICS_EXE}" \
      --input "\$IN_PATH" \
      --raw_layer "${raw_layer}" \
      --normalized_layer "${normalized_layer}" \
      --max_cells "${max_cells}" \
      --random_state "${random_state}" \
      --output "\$OUT_PATH"
    """
}

process compute_segmentation_quality_metrics {
    tag "metrics_segmentation_quality"
    publishDir PUBLISH_DIR, mode: 'copy'

    input:
    path adata_in

    output:
    path "metrics_segmentation_quality.h5ad", emit: metrics

    script:
    def layer = params.metric_marker_layer ?: 'counts'
    def nuclear_markers = params.nuclear_markers ?: 'null'
    def exclusive_pairs = params.exclusive_pairs ?: 'null'
    def n_neighbors = params.containsKey('n_neighbors') ? params.n_neighbors : 6
    def max_cells = params.containsKey('metric_max_cells') ? params.metric_max_cells : 5000
    def random_state = params.containsKey('metric_random_state') ? params.metric_random_state : 0
    """
    set -euo pipefail
    export VIASH_ENGINE_TYPE=${VIASH_ENGINE}

    "${METRICS_SEGMENTATION_QUALITY_EXE}" \
      --input "\$(realpath '${adata_in}')" \
      --layer "${layer}" \
      --nuclear_markers "${nuclear_markers}" \
      --exclusive_pairs "${exclusive_pairs}" \
      --n_neighbors "${n_neighbors}" \
      --max_cells "${max_cells}" \
      --random_state "${random_state}" \
      --output "\$PWD/metrics_segmentation_quality.h5ad"
    """
}

process compute_technical_covariate_metrics {
    tag "metrics_technical_covariates"
    publishDir PUBLISH_DIR, mode: 'copy'

    input:
    path adata_in

    output:
    path "metrics_technical_covariates.h5ad", emit: metrics

    script:
    def layer = params.metric_normalized_layer ?: 'normalized'
    def covariate = params.metric_covariate ?: 'area'
    def n_pcs = params.containsKey('metric_n_pcs') ? params.metric_n_pcs : 10
    def max_cells = params.containsKey('metric_max_cells') ? params.metric_max_cells : 5000
    def random_state = params.containsKey('metric_random_state') ? params.metric_random_state : 0
    """
    set -euo pipefail
    export VIASH_ENGINE_TYPE=${VIASH_ENGINE}

    "${METRICS_TECHNICAL_COVARIATES_EXE}" \
      --input "\$(realpath '${adata_in}')" \
      --layer "${layer}" \
      --covariate "${covariate}" \
      --n_pcs "${n_pcs}" \
      --max_cells "${max_cells}" \
      --random_state "${random_state}" \
      --output "\$PWD/metrics_technical_covariates.h5ad"
    """
}

process compute_signal_background_metrics {
    tag "metrics_signal_background"
    publishDir PUBLISH_DIR, mode: 'copy'

    input:
    path corrected_image
    path segmentation_labels
    path corrected_markers
    path adata_qc

    output:
    path "metrics_signal_background.h5ad", emit: metrics

    script:
    def dilation_iterations = params.containsKey('dilation_iterations') ? params.dilation_iterations : 3
    def norm_percentile = params.containsKey('norm_percentile') ? params.norm_percentile : 95.0
    def random_state = params.containsKey('metric_random_state') ? params.metric_random_state : 0
    """
    set -euo pipefail
    export VIASH_ENGINE_TYPE=${VIASH_ENGINE}

    "${METRICS_SIGNAL_BACKGROUND_EXE}" \
      --input_image "\$(realpath '${corrected_image}')" \
      --input_labels "\$(realpath '${segmentation_labels}')" \
      --input_markers "\$(realpath '${corrected_markers}')" \
      --input_qc "\$(realpath '${adata_qc}')" \
      --dilation_iterations "${dilation_iterations}" \
      --norm_percentile "${norm_percentile}" \
      --random_state "${random_state}" \
      --output "\$PWD/metrics_signal_background.h5ad"
    """
}

workflow run_wf {
    input_zarr = params.zarr ?: params.input_zarr

    if (input_zarr == null) {
        error "No input Zarr provided. Use --zarr or --input_zarr."
    }

    zarr_ch = Channel.value(file(input_zarr))

    corrected = image_correction(zarr_ch)
    segmentation = segment_cells(corrected.image, zarr_ch)
    aggregated = aggregate_proteins(corrected.image, segmentation.labels, corrected.markers)
    qc = qc_filter(aggregated.adata)
    norm = normalize(qc.adata)
    basic_metrics = compute_basic_metrics(norm.adata)
    marker_metrics = compute_marker_distribution_metrics(norm.adata)
    normalization_metrics = compute_normalization_diagnostics(norm.adata)
    segmentation_quality_metrics = compute_segmentation_quality_metrics(norm.adata)
    technical_covariate_metrics = compute_technical_covariate_metrics(norm.adata)
    signal_background_metrics = compute_signal_background_metrics(
        corrected.image, segmentation.labels, corrected.markers, qc.adata
    )

    basic_metrics.metrics.view { "Done. Basic metrics written to: ${it}" }
    marker_metrics.metrics.view { "Done. Marker distribution metrics written to: ${it}" }
    normalization_metrics.metrics.view { "Done. Normalization diagnostics written to: ${it}" }
    segmentation_quality_metrics.metrics.view { "Done. Segmentation quality metrics written to: ${it}" }
    technical_covariate_metrics.metrics.view { "Done. Technical covariate metrics written to: ${it}" }
    signal_background_metrics.metrics.view { "Done. Signal/background metrics written to: ${it}" }
}
