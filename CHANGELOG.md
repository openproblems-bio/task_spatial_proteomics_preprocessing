# task_spatial_proteomics_preprocessing

## Unreleased

Initial port of the spatial proteomics preprocessing benchmark onto the
OpenProblems task structure. Developed in
`juliusstein77/task_isp_preprocessing`; only the proteomics surface was
carried across.

### Components

* Five method namespaces, each with interchangeable alternatives sharing one
  argument contract: `methods_image_correction`, `methods_proteomics_segmentation`,
  `methods_protein_aggregation`, `methods_proteomics_qc`,
  `methods_proteomics_normalization`.
* Six metric families in `metrics_proteomics`, wrapping
  [`sprotmetrics`](https://github.com/theislab/sprotmetrics) `0.1.0`.

### Segmentation methods ported from task_ist_preprocessing

* Five segmentation components adapted from
  `openproblems-bio/task_ist_preprocessing`: `cellpose_segmentation`
  (Cellpose 3: cyto, cyto2, cyto3 or nuclei), `cellposev4_segmentation`
  (Cellpose-SAM), `stardist_segmentation`, `watershed_txsim_segmentation`
  and `binning_segmentation`. The algorithms and default hyperparameters are
  kept. I/O is rewritten to this task's `.npy` contract (gapless int32
  labels and a metadata JSON). The txsim watershed and binning logic is
  inlined, so these components no longer depend on txsim, squidpy or
  rasterio. Integer-valued images (stored as float32 by image correction)
  are cast back to uint16 and take txsim's exact code path. Against
  `txsim.preprocessing.segment_watershed` with scikit-image 0.26 the masks
  were pixel-identical on 6 synthetic images × 6 option sets (IST defaults,
  equalize_hist, median blur, global Otsu, bg factor 0.9, compact
  watershed). Genuinely float images are scaled to [0, 1] first, which can
  shift a few boundary pixels. `custom_segmentation` was not ported because
  `reference_segmentation` already covers it.
* The deep-learning components carry the same GPU directive labels as in
  IST (`gpuhighmem` for both Cellpose components, `gpu` for StarDist). As in
  IST, these labels are not defined in `_viash.yaml`'s `config_mods`, so the
  Nextflow profile has to map them to GPU resources. StarDist and Cellpose
  weights are downloaded into the images at build time.
* `remove_small_objects` / `remove_small_holes` calls go through a
  version-agnostic helper. It uses `max_size` where scikit-image provides it
  (0.26+, where `min_size` / `area_threshold` are deprecated), with the same
  numeric value, which gives the same masks without the warning. Outputs of
  `basic_threshold_segmentation` and `watershed_expand_segmentation` were
  checked to be unchanged.
* Segmentation methods can select channels by marker name.
  `--input_markers` (the image-correction sidecar) together with
  `--nuclear_markers` and `--membrane_markers` replaces the IST components'
  hard-coded `image[0]`. Several names are percentile-scaled and
  max-projected. `--nuclei_channel` remains the fallback. Both Cellpose
  components run in two-channel [membrane, nuclear] mode when membrane
  markers are given.
* The segmentation API now carries only the arguments every method
  understands. The threshold and watershed arguments moved into
  `basic_threshold_segmentation` and `watershed_expand_segmentation`, whose
  behaviour is unchanged.
* `run_benchmark` runs all seven image-derived segmenters by default. It
  forwards the marker sidecar to segmentation and takes new
  `--segmentation_nuclear_markers` and `--membrane_markers` arguments.

### Conformance to the OpenProblems structure

* `src/api/file_score.yaml` describes what a metric component actually writes:
  `metric_ids`, `metric_values` and the metric family. It no longer asks for
  `dataset_id`, `normalization_id` and `method_id`, which a metric here cannot
  know - the result is produced by a chain of five methods and the file a metric
  reads carries no record of the four steps before the last. Provenance is
  attached by the workflow instead, which is the only place that has it. The six
  metric API contracts merge this file into their `--output`, so the contract is
  wired rather than documentary.
* Every metric component declares the `info.metrics` block that
  `common/schemas/task_metric.yaml` requires - 66 metrics across six families,
  generated from the `sprotmetrics` registry by
  `scripts/create_component/generate_metric_info.py` rather than maintained by
  hand. Metrics with no direction are deliberately omitted from that block: the
  schema requires `maximize`, and has no way to express "descriptive", so
  declaring them would hand the scoring stage a ranking criterion they cannot
  support. They are still computed and still appear in the score file.
* The workflow moves from `src/workflows/run_benchmark_isp/main_isp.nf` to the
  conventional `src/workflows/run_benchmark/main.nf`.

### Workflow

* `run_benchmark` is now built on the OpenProblems `runEach` / `fromState` /
  `toState` idiom, replacing the hand-written processes that shelled out to
  `target/executable/...` with a `case` block per stage. A state threads down
  the chain carrying the files produced so far and a `steps` list recording
  which component produced each, and the workflow ends in a scoring stage.

  Two things had to change in the components first, and neither was visible
  while the workflow called executables:

  - No component declared a `nextflow` runner, only `executable`. Viash cannot
    emit a component as a Nextflow module without it, so composing modules was
    not an option at all. The failure does not say so: `viash ns build` reports
    the workflow's 21 dependencies as *missing*.
  - No output file argument declared an `example`. Viash takes the extension of
    an auto-generated output filename from `example`, so without one the file
    gets no extension - and `np.save` appends `.npy` to a path that lacks it.
    The component wrote `...output_image.npy` while the module looked for
    `...output_image` and reported it missing. The old workflow passed explicit
    paths that already carried the extension, which is why this never showed.

* One invocation now runs many method chains rather than one, and the choice of
  which is an argument rather than an outer bash loop. `--default_methods` names
  the baseline, one method per step, and a method that is not a default is only
  combined with the defaults - so each non-baseline chain differs from the
  baseline in exactly one step, which is what makes a difference in a metric
  attributable to a single choice. With the shipped defaults that is 10 chains
  rather than the 162 a cross-product would give.

* Per-method parameters are no longer workflow arguments. The 40-argument
  config collapses to the method-selection lists plus the three metric
  arguments that depend on the antibody panel; everything else comes from the
  component defaults, or from `--method_parameters_yaml`, which can also sweep
  one parameter across values. Sweeps expand only against the baseline, for the
  same reason non-default methods do.

* `--viash_engine` is gone. It existed because every process hardcoded
  `export VIASH_ENGINE_TYPE=docker`; the workflow now composes Viash modules, so
  the engine is a Nextflow profile (`-profile docker`, `-profile singularity`,
  or none to run against the local environment).

* The run publishes `score_uns.yaml`, `method_configs.yaml`,
  `metric_configs.yaml`, `dataset_uns.yaml` and `task_info.yaml`. Scores carry
  the dataset, the full chain, the normalization, the metric component and the
  run id of every step - the provenance the hand-written workflow could not
  produce, and the reason `file_score.yaml` no longer asks metric components
  for it.

  Verified end to end on a 512x512 crop of exemplar-001, under Viash's native
  engine: 10 chains x 6 metric families = 60 score records over 82 distinct
  metric ids, each chain differing from the baseline in exactly one step. The
  axes separate as they do on the full image - `n_cells` 239 for threshold
  against 293 for watershed, `qc_pass_rate` 1.00 / 0.91 / 0.55 across the three
  QC methods, and `variance_explained_by_covariate_counts` 0.788 for `sum`
  against 0.035 for `mean`. Chains differing only in normalization report
  identical segmentation metrics, as they must when normalization is per-cell.

* `src/api/file_spatial_proteomics_dataset.yaml` gives the input a contract.
  The workflow takes a SpatialData zarr and the repository described that only
  in prose.

* `scripts/run_benchmark/*.sh` run this task rather than the template's
  train/test/solution split, which was left over from the OpenProblems task
  template and referred to files this task has never produced.

### Known gaps

* `dataset_uns.yaml` carries little more than the run id. The input is a
  SpatialData zarr rather than an h5ad, so `utils/extract_uns_metadata` cannot
  read metadata off it, and there is no datasets workflow that would write that
  metadata into the store in the first place. Real dataset metadata belongs with
  the `process_datasets` work below.
* Test resources are not yet on S3. `_viash.yaml` carries no `test_resources`
  entry, because naming a bucket path that does not exist makes CI's
  resource-sync step fail and skip every component test - red while testing
  nothing. The intended entry is recorded there as a comment. `scripts/create_resources/isp/` builds the input locally from the
  public MCMICRO download in the meantime.
* There is no `workflows/process_datasets`. The template's version implements a
  train/test/solution split for label projection, which is not the shape of this
  task - the input is a SpatialData store, and preprocessing is what is being
  benchmarked rather than prediction. It was removed rather than left dangling
  against deleted API files. A real dataset processor belongs with the S3
  resources work above.
