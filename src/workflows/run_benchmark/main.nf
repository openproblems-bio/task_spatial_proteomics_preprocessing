include { checkItemAllowed } from "${meta.resources_dir}/helper.nf"

/*
 * Spatial proteomics preprocessing benchmark.
 *
 *   SpatialData zarr
 *     -> image correction
 *     -> segmentation
 *     -> cell-by-protein aggregation
 *     -> quality control
 *     -> normalization
 *     -> metrics
 *
 * Each step is a `runEach` over the components registered for it. A state
 * travels down the chain carrying the files produced so far and a `steps` list
 * recording which component produced each one, so the score file can say which
 * chain a metric value came from.
 *
 * Which chains are run is decided by `expandChannelWithParameterSets`, below:
 * `--default_methods` names one method per step, and a method that is not a
 * default is only ever combined with the defaults. The run therefore varies one
 * step at a time rather than taking the cross-product, which is what makes a
 * difference in a metric attributable to a single choice.
 */

// Step types that name a preprocessing choice, as opposed to the dataset and
// the metric. Used to build the method id from a chain.
METHOD_STEP_TYPES = ["image_correction", "segmentation", "aggregation", "qc", "normalization"]

workflow auto {
  findStates(params, meta.config)
    | meta.workflow.run(
      auto: [publish: "state"]
    )
}

workflow run_wf {
  take:
  input_ch

  main:

  /****************************************
   *      INITIALIZE DATA STRUCTURES      *
   ****************************************/
  init_ch = input_ch
    | map { id, state ->
      // optional per-method parameter sets and sweeps
      def method_parameters = null
      if (state.method_parameters_yaml) {
        def yaml_blob = readYaml(state.method_parameters_yaml)
        if (yaml_blob instanceof Map && yaml_blob.containsKey('parameters')) {
          method_parameters = yaml_blob.parameters
        } else if (yaml_blob instanceof Map) {
          method_parameters = yaml_blob
        }
      }

      def new_state = state + [
        orig_id: id,
        steps: [[type: "dataset", dataset_id: id]],
        method_parameters: method_parameters
      ]
      [id, new_state]
    }

  /****************************************
   *        RUN IMAGE CORRECTION          *
   ****************************************/
  corr_methods = [
    no_image_correction,
    basic_image_correction,
    flatfield_image_correction
  ]

  corr_ch = init_ch
    | expandChannelWithParameterSets(corr_methods, "corr", "image_correction_methods")
    | runEach(
      components: corr_methods,
      filter: { id, state, comp ->
        comp.config.name == state.current_method_id
      },
      fromState: { id, state, comp ->
        def args = [input_zarr: state.input_zarr]
        if (state.image_key) {
          args.image_key = state.image_key
        }
        args + state.current_method_args
      },
      toState: { id, out_dict, state, comp ->
        recordStep(state, "image_correction", id) + [
          output_image: out_dict.output_image,
          output_markers: out_dict.output_markers
        ]
      }
    )

  /****************************************
   *          RUN SEGMENTATION            *
   ****************************************/
  segm_methods = [
    basic_threshold_segmentation,
    watershed_expand_segmentation,
    watershed_txsim_segmentation,
    binning_segmentation,
    cellpose_segmentation,
    cellposev4_segmentation,
    stardist_segmentation,
    reference_segmentation
  ]

  segm_ch = corr_ch
    | expandChannelWithParameterSets(segm_methods, "segm", "segmentation_methods")
    | runEach(
      components: segm_methods,
      filter: { id, state, comp ->
        comp.config.name == state.current_method_id
      },
      fromState: { id, state, comp ->
        def args = [:]
        if (comp.config.name == "reference_segmentation") {
          // reads the dataset's own labels rather than the corrected image
          args.input_zarr = state.input_zarr
          if (state.labels_key) {
            args.labels_key = state.labels_key
          }
        } else {
          args.input_image = state.output_image
          // marker names let methods pick channels by name, not index
          args.input_markers = state.output_markers
          if (state.segmentation_nuclear_markers) {
            args.nuclear_markers = state.segmentation_nuclear_markers
          }
          if (state.membrane_markers) {
            args.membrane_markers = state.membrane_markers
          }
        }
        args + state.current_method_args
      },
      toState: { id, out_dict, state, comp ->
        recordStep(state, "segmentation", id) + [
          output_labels: out_dict.output_labels,
          output_segmentation_metadata: out_dict.output_metadata
        ]
      }
    )

  /****************************************
   *       RUN PROTEIN AGGREGATION        *
   ****************************************/
  aggr_methods = [
    mean_intensity_aggregation,
    sum_intensity_aggregation,
    median_intensity_aggregation
  ]

  aggr_ch = segm_ch
    | expandChannelWithParameterSets(aggr_methods, "aggr", "aggregation_methods")
    | runEach(
      components: aggr_methods,
      filter: { id, state, comp ->
        comp.config.name == state.current_method_id
      },
      fromState: { id, state, comp ->
        [
          input_image: state.output_image,
          input_labels: state.output_labels,
          input_markers: state.output_markers
        ] + state.current_method_args
      },
      toState: { id, out_dict, state, comp ->
        recordStep(state, "aggregation", id) + [
          output_aggregation: out_dict.output
        ]
      }
    )

  /****************************************
   *            RUN QC FILTER             *
   ****************************************/
  qc_methods = [
    no_qc_filter,
    basic_qc_filter,
    adaptive_quantile_qc_filter
  ]

  qc_ch = aggr_ch
    | expandChannelWithParameterSets(qc_methods, "qc", "qc_methods")
    | runEach(
      components: qc_methods,
      filter: { id, state, comp ->
        comp.config.name == state.current_method_id
      },
      fromState: { id, state, comp ->
        [input: state.output_aggregation] + state.current_method_args
      },
      toState: { id, out_dict, state, comp ->
        recordStep(state, "qc", id) + [
          output_qc: out_dict.output
        ]
      }
    )

  /****************************************
   *          RUN NORMALIZATION           *
   ****************************************/
  norm_methods = [
    no_normalization,
    total_log1p_normalization,
    arcsinh_normalization
  ]

  norm_ch = qc_ch
    | expandChannelWithParameterSets(norm_methods, "norm", "normalization_methods")
    | runEach(
      components: norm_methods,
      filter: { id, state, comp ->
        comp.config.name == state.current_method_id
      },
      fromState: { id, state, comp ->
        [input: state.output_qc] + state.current_method_args
      },
      toState: { id, out_dict, state, comp ->
        recordStep(state, "normalization", id) + [
          output_normalization: out_dict.output
        ]
      }
    )

  /****************************************
   *                METRICS               *
   ****************************************/
  metrics = [
    basic,
    marker_distribution,
    normalization_diagnostics,
    segmentation_quality,
    technical_covariates,
    signal_background
  ]

  metric_ch = norm_ch
    | runEach(
      components: metrics,
      id: { id, state, comp ->
        id + "/metric_" + comp.config.name
      },
      fromState: { id, state, comp ->
        // signal_background is the only metric computed from pixels rather
        // than from the cell table, so it takes the image, the masks and the
        // QC annotation instead of the normalized matrix.
        if (comp.config.name == "signal_background") {
          return [
            input_image: state.output_image,
            input_labels: state.output_labels,
            input_markers: state.output_markers,
            input_qc: state.output_qc,
            dilation_iterations: state.dilation_iterations
          ]
        }

        def args = [input: state.output_normalization]
        if (comp.config.name == "segmentation_quality") {
          // panel-specific; the metrics stay NaN rather than guess when unset
          if (state.nuclear_markers) {
            args.nuclear_markers = state.nuclear_markers
          }
          if (state.exclusive_pairs) {
            args.exclusive_pairs = state.exclusive_pairs
          }
        }
        return args
      },
      toState: { id, out_dict, state, comp ->
        state + [
          steps: state.steps + [[
            type: "metric",
            component_id: comp.config.name,
            component_variant: comp.config.name,
            run_id: id
          ]],
          output_metric: out_dict.output
        ]
      }
    )

    // read uns['metric_ids'] and uns['metric_values'] back out of the h5ad
    | extract_uns_metadata.run(
      key: "extract_uns_scores",
      fromState: [input: "output_metric"],
      args: [
        uns_length_cutoff: 100
      ],
      toState: { id, output, state ->
        state + [
          score_uns: readYaml(output.output).uns
        ]
      }
    )

    | joinStates { ids, states ->
      def score_uns = states.collect { state ->
        def method_steps = state.steps.findAll { it.type in METHOD_STEP_TYPES }
        def metric_step = state.steps.find { it.type == "metric" }
        def norm_step = method_steps.find { it.type == "normalization" }

        [
          dataset_id: state.orig_id,
          // the whole chain, since no single component owns the result
          method_id: method_steps.collect { it.component_variant }.join("+"),
          method_ids: method_steps.collect { it.component_id },
          normalization_id: norm_step?.component_variant,
          metric_component: metric_step?.component_id,
          steps: state.steps,
          metric_ids: state.score_uns['metric_ids'],
          metric_values: state.score_uns['metric_values']
        ]
      }

      def score_uns_file = tempFile("score_uns.yaml")
      score_uns_file.write(toYamlBlob(score_uns))

      ["output", [output_scores: score_uns_file]]
    }

  /******************************
   * GENERATE OUTPUT YAML FILES *
   ******************************/
  meta_ch = init_ch
    | map { id, state ->
      [id, state + ["_meta": [join_id: id]]]
    }
    | joinStates { ids, states ->
      // The input is a SpatialData zarr, not an h5ad, so there is no `uns` for
      // utils/extract_uns_metadata to read and the dataset metadata has to be
      // assembled from the run itself. Anything richer than this belongs in a
      // datasets workflow that writes it into the store.
      def dataset_uns = states.collect { state ->
        [
          dataset_id: state.orig_id,
          dataset_name: state.orig_id,
          image_key: state.image_key,
          labels_key: state.labels_key
        ]
      }
      def dataset_uns_file = tempFile("dataset_uns.yaml")
      dataset_uns_file.write(toYamlBlob(dataset_uns))

      def methods = corr_methods + segm_methods + aggr_methods + qc_methods + norm_methods
      def method_configs_file = tempFile("method_configs.yaml")
      method_configs_file.write(toYamlBlob(methods.collect { it.config }))

      def metric_configs_file = tempFile("metric_configs.yaml")
      metric_configs_file.write(toYamlBlob(metrics.collect { it.config }))

      def viash_file = meta.resources_dir.resolve("_viash.yaml")

      def new_state = [
        output_dataset_info: dataset_uns_file,
        output_method_configs: method_configs_file,
        output_metric_configs: metric_configs_file,
        output_task_info: viash_file,
        _meta: states[0]._meta
      ]

      ["output", new_state]
    }

  output_ch = metric_ch
    | mix(meta_ch)
    | joinStates { ids, states ->
      def mergedStates = states.inject([:]) { acc, m -> acc + m }
      [ids[0], mergedStates]
    }

  emit:
  output_ch
}

/**
 * Remove specified keys from a map.
 */
def removeKeys(map, keys) {
  def result = map.clone()
  keys.each { key -> result.remove(key) }
  return result
}

/**
 * Append one preprocessing step to the state's provenance and clear the
 * per-run method selection, so the next step starts from a clean state.
 */
def recordStep(state, step_type, run_id) {
  removeKeys(state, ["current_method_id", "current_method_variant", "current_method_args"]) + [
    steps: state.steps + [[
      type: step_type,
      component_id: state.current_method_id,
      component_variant: state.current_method_variant,
      run_id: run_id
    ]]
  ]
}

/**
 * Expand a channel into one copy per method, and per parameter-set variant.
 *
 * Two things limit the expansion, both so that a chain differs from the
 * baseline in at most one place:
 *
 *   - a method that is not in `--default_methods` is skipped once the state
 *     already carries a non-default choice;
 *   - a parameter sweep is expanded only while the state is still the baseline.
 *
 * @param methods List of components for this step
 * @param method_type Short tag used in the generated run ids
 * @param selected_methods_key State key holding the methods selected for this step
 * @return Workflow that takes a channel and emits expanded [id, state] tuples
 */
def expandChannelWithParameterSets(methods, method_type, selected_methods_key) {
  workflow expand_ch_wf {
    take: input_ch

    main:
    output_ch = input_ch
      | flatMap { id, state ->
        def selected_methods = state[selected_methods_key]
        def default_methods = state.default_methods ?: []
        def results = []

        methods.each { comp ->
          def comp_name = comp.config.name

          if (!checkItemAllowed(comp_name, selected_methods, null, selected_methods_key, "NA")) {
            return
          }

          // does the state already differ from the baseline?
          def has_non_default = state.steps.any { step ->
            step.component_id != null &&
            (!(step.component_id in default_methods) ||
             (step.component_variant != null && step.component_id != step.component_variant))
          }

          def is_non_default_method = !(comp_name in default_methods)
          if (has_non_default && is_non_default_method) {
            return
          }

          def spec = state.method_parameters ? state.method_parameters[comp_name] : null
          def default_args = spec ? (spec['default'] ?: [:]) : [:]

          results << [
            "${id}/${method_type}_${comp_name}".toString(),
            state + [
              current_method_id: comp_name,
              current_method_variant: comp_name,
              current_method_args: default_args
            ]
          ]

          def sweep = spec ? spec['sweep'] : null
          if (sweep instanceof Map && !has_non_default) {
            sweep.each { parName, parVals ->
              if (!(parVals instanceof Collection)) parVals = [parVals]
              parVals.each { val ->
                def variant_name = "${comp_name}_${parName}_${val}".toString()
                results << [
                  "${id}/${method_type}_${variant_name}".toString(),
                  state + [
                    current_method_id: comp_name,
                    current_method_variant: variant_name,
                    current_method_args: default_args + [(parName): val]
                  ]
                ]
              }
            }
          }
        }

        return results
      }

    emit:
    output_ch
  }
  return expand_ch_wf.cloneWithName("expand_channel_for_" + method_type)
}
