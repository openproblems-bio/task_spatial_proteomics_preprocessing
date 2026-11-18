# Developer guide

How to work on the spatial proteomics (ISP) preprocessing benchmark: the layout,
how to add a method, and how to run the tests.

> This repository began as a fork of the OpenProblems iST preprocessing task.
> A good deal of the original imaging-based spatial *transcriptomics* code is
> still present under `src/` (`methods_transcript_assignment/`,
> `methods_segmentation/`, `methods_cell_type_annotation/`,
> `methods_expression_correction/`, `metrics/`, `src/workflows/run_benchmark/`).
> None of it is used by the ISP workflow. When in doubt, the ISP surface is
> exactly what `src/workflows/run_benchmark/config.vsh.yaml` lists under
> `dependencies`.

## Tooling

| Tool | Version | Purpose |
|---|---|---|
| [Viash](https://viash.io) | 0.9.4 (pinned in `_viash.yaml`) | Turns `config.vsh.yaml` + `script.py` into a standalone executable |
| Nextflow | tested with 24.04.3 | Runs the workflow |
| Docker | any recent | Default execution engine for components |

Every component also declares a `native` engine, so the suite can be built and
tested against a local Python environment when Docker is unavailable. See
[Running the tests](#running-the-tests).

## Layout

```text
src/
  api/                    Shared component contracts (arguments + test wiring)
  methods_image_correction/
  methods_proteomics_segmentation/
  methods_protein_aggregation/
  methods_proteomics_qc/
  methods_proteomics_normalization/
  metrics_proteomics/
  testing/                Synthetic fixtures and the shared test scripts
  workflows/run_benchmark_isp/
                          main.nf, config.vsh.yaml, test.sh
scripts/thesis_benchmark/ Drivers for the 15 published runs, plus figure/table generation
docs/metrics_catalogue.csv What every metric measures, its direction, and its origin
```

A component is a directory with two files:

- `config.vsh.yaml` — merges an API contract via `__merge__`, then adds the
  component's name, description, and engine dependencies.
- `script.py` — reads a `par` dict that Viash injects, writes the declared
  outputs.

The API contract in `src/api/comp_*.yaml` owns the argument list, so every
method in a namespace accepts the same flags. That is what lets the workflow
swap one method for another without changing its call.

## The ISP pipeline

```text
SpatialData .zarr
  -> image correction   corrected_image.npy (C x Y x X float32) + corrected_markers.json
  -> segmentation       segmentation_labels.npy (Y x X int32) + segmentation_metadata.json
  -> aggregation        cell_by_protein.h5ad
  -> QC                 cell_by_protein_qc.h5ad
  -> normalization      cell_by_protein_norm.h5ad
  -> metrics            metrics_*.h5ad
```

Two contracts are worth knowing before you change anything:

- **Segmentation labels are a gapless `1..n` integer range**, so
  `max(labels) == n_objects`. Image-derived methods must renumber after any
  step that drops objects (`relabel_sequential`). Reference segmentation is the
  deliberate exception: dataset-provided label ids are passed through untouched
  so they stay comparable to the dataset's own cell ids.
- **QC annotates, it does not filter.** QC components add `obs['passed_QC']`,
  `obs['qc_filter_reason']` and `uns['qc']`, and keep every cell. Normalization
  runs on all of them, which is harmless because every normalization method is
  per-cell: which cells are present cannot change the values of the cells that
  remain. Metrics are the only population-level statistics in the pipeline, so
  they are where QC has to be accounted for, and each metric component reports
  its values over all cells and again over the QC-passing subset with a
  `_qc_pass` suffix. A new metric component should do the same, or the QC module
  will look inert through it.

## Adding a method

Say you are adding a segmentation method.

1. **Create the component.**

   ```bash
   mkdir -p src/methods_proteomics_segmentation/my_method
   ```

   `config.vsh.yaml`:

   ```yaml
   __merge__: /src/api/comp_method_proteomics_segmentation.yaml

   name: my_method
   label: "My method"
   summary: "One line."
   description: |
     What it does and what it is a baseline for.

   resources:
     - type: python_script
       path: script.py

   engines:
     - type: docker
       image: openproblems/base_python:1
       setup:
         - type: python
           pypi:
             - numpy
             - scikit-image
     - type: native

   runners:
     - type: executable
   ```

   Keep the `native` engine: without it the component is invisible to
   `viash ns test --engine native` and silently drops out of the test suite.

2. **Write `script.py`.** Read from `par`, write the declared outputs, and end
   with the `## VIASH START` / `## VIASH END` block giving defaults for direct
   execution. Prefer raising on an unrecognised option over falling back to a
   default — a silently ignored argument is reported as a successful run.

3. **Add any new arguments to the API contract**, not just to your component.
   The workflow calls every method in a namespace with the same flags, so an
   argument that only one component declares will be dropped for the others.
   Viash 0.9 does *not* reject unknown arguments; it warns and treats them as
   positional, so this failure mode is silent.

4. **Register it with the workflow** in three places:
   - `src/workflows/run_benchmark/config.vsh.yaml` under `dependencies`
   - `main.nf`: an `*_EXE` constant and a branch in the stage's `case`
   - `src/workflows/run_benchmark/test.sh`: the build list

5. **Run the tests** (below). The shared test script for the namespace picks up
   the new component automatically; extend it if the method needs a new
   assertion.

## Running the tests

Each API contract wires in a shared test script from `src/testing/`, so every
component in a namespace is covered by one file:

| Namespace | Test |
|---|---|
| `methods_image_correction` | `test_image_correction.py` |
| `methods_proteomics_segmentation` | `test_segmentation.py` |
| `methods_protein_aggregation` (arrays) | `test_aggregation.py` |
| `methods_protein_aggregation` (zarr) | `test_aggregation_zarr.py` |
| `methods_proteomics_qc` | `test_qc.py` |
| `methods_proteomics_normalization` | `test_normalization.py` |
| `metrics_proteomics` (matrix-based) | `test_metrics.py` |
| `metrics_proteomics/signal_background` (pixel-based) | `test_signal_background.py` |

The tests build their own inputs with `src/testing/isp_fixtures.py` — a 64x64,
3-channel synthetic image, a matching label mask, and small AnnData objects.
Nothing is downloaded, so the whole suite runs in well under a minute and needs
no access to the MCMICRO dataset.

With Docker:

```bash
viash ns test --query "methods_proteomics_segmentation"
```

Without Docker, against a local environment that provides the packages the
components declare (`numpy`, `scipy`, `scikit-image`, `pandas`, `anndata`,
`spatialdata`, `zarr`):

```bash
viash ns test --engine native --query "methods_proteomics_segmentation"
```

Everything at once:

```bash
for ns in methods_image_correction methods_proteomics_segmentation \
          methods_protein_aggregation methods_proteomics_qc \
          methods_proteomics_normalization metrics_proteomics; do
  viash ns test --query "$ns"
done
```

Note that `viash ns test` exits 0 even when individual tests fail, so read the
`result` column rather than the exit status. A namespace that reports no rows
at all usually means no component in it declares the engine you asked for.

## Running the workflow

```bash
ZARR_PATH=/path/to/exemplar-001.zarr \
bash src/workflows/run_benchmark/test.sh
```

`test.sh` rebuilds every component and then runs `main.nf`. It must be
launched from the repository root: `main.nf` resolves component executables
relative to `workflow.launchDir`.

`VIASH_ENGINE` (default `docker`) selects the engine both for the build and for
the run; the two must agree, which is why one variable drives both. Set it to
`native` to run without Docker:

```bash
VIASH_ENGINE=native ZARR_PATH=/path/to/exemplar-001.zarr \
bash src/workflows/run_benchmark/test.sh
```

The dataset itself is not in the repository and the public download is not in
the form the workflow wants; `scripts/create_resources/isp/prepare_exemplar_001.sh`
builds it. See **Data availability** in [README.md](README.md).

Method choice and parameters come from environment variables — see the table in
[README.md](README.md). To vary one module at a time across the published
comparison set, use `scripts/thesis_benchmark/run_15_benchmark_runs.sh`.

## Known rough edges

Worth knowing before you plan larger changes:

- `main.nf` is hand-written Nextflow with a `case` block per stage, rather
  than the OpenProblems `runEach` / `fromState` idiom used by the legacy iST
  workflow in `src/workflows/run_benchmark/main.nf`. One run covers exactly one
  method combination; sweeps are done by an outer bash loop.
- There is no scoring stage in the workflow. Result tables are assembled
  afterwards by `scripts/thesis_benchmark/collect_benchmark_results.py`.
- Image components load the whole image with `np.asarray(...)`, with no dask or
  multiscale handling, and infer channel order from a shape heuristic. The
  corrected image is written as a dense float32 `.npy`, which for exemplar-001
  is 371 MB per run.
- Several components are superseded but retained: `normalize_total_intensity`,
  `basic_protein_aggregation`, `basic_protein_aggregation_arrays`.
  `metrics_proteomics/cluster_separability` and the two
  `methods_proteomics_cell_type_annotation` components are implemented but not
  wired into the workflow.
