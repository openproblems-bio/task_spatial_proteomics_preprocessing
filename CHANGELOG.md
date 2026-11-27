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

### Conformance to the OpenProblems structure

* `src/api/file_score.yaml` now requires `dataset_id`, `normalization_id` and
  `method_id` in `uns`. These were commented out in the source repository, so
  scores carried no provenance and could not be attributed to a run.
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

### Known gaps

* `run_benchmark/main.nf` is still the hand-written form inherited from the
  source repository - one process per stage with a `case` block - rather than
  the `runEach` / `fromState` idiom the template uses. Porting it is what will
  supply the `method_id` and `dataset_id` that `file_score.yaml` now requires.
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
