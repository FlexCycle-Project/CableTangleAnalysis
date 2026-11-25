# CableTangleAnalysis

Julia analysis of CableSim tangle metrics.  The package produces the paper
figures and tables:

- knot-type distribution (bar chart),
- PCA projection of the five geometric metrics,
- correlation heatmap with FDR-corrected p-values,
- cut-strategy comparisons for threading reduction and cuts needed,
- a PCA-derived complexity score and its reduction after cutting,
- booktabs LaTeX tables for the PCA and complexity results.

`STATISTICS.md` documents the statistical methods (paired Wilcoxon tests,
FDR correction, effect sizes, confidence intervals, symmetry checks).

## Prerequisites

Julia >= 1.6.  Install dependencies once:

```bash
julia --project=cabletangleanalysis -e 'using Pkg; Pkg.instantiate()'
```

## Run the dataset analysis

The entry point for a CableSim dataset is `scripts/run_dataset.jl`:

```bash
julia --project=cabletangleanalysis \
  cabletangleanalysis/scripts/run_dataset.jl <dataset-root> [variant ...]
```

A variant may be a directory name below the dataset root (`zero`), a path
(`results/run/zero`), or `label=path` when the output label should differ from
the directory name.  Datasets whose tables sit directly in the root, such as
`results/dataset_final`, use the label form for a single variant:

```bash
julia --project=cabletangleanalysis \
  cabletangleanalysis/scripts/run_dataset.jl \
  results/dataset_final single=results/dataset_final
```

Without a variant list the driver runs `zero`, `current`, and `high`.

Each variant directory must contain:

```text
catalog.tsv
cut_eval.tsv
cuts_to_zero.tsv
cut_eval_stage_metrics.tsv
cuts_to_zero_stage_metrics.tsv
```

Outputs are written to `<dataset-root>/analysis_julia/<variant>/`:

```text
figures/knot_distribution.pdf
figures/pca/pca2d_combined.pdf
figures/correlation_heatmap.pdf
figures/cut_eval/cuts_histogram.pdf
figures/cut_eval/cut_eval_comparison.pdf
figures/cut_eval/threading_reduction_first_cut.pdf  (+ _paired.pdf)
figures/cut_eval/threading_reduction_resolved.pdf   (+ _paired.pdf)
figures/cut_eval/cuts_needed_all.pdf                (+ _paired.pdf)
figures/cut_eval/cuts_needed_resolved.pdf           (+ _paired.pdf)
figures/cut_eval/cut_eval_comparison_complexity.pdf
figures/cut_eval/complexity_reduction_resolved.pdf
tables/pca_results.tsv
tables/pca_scores.tsv
tables/*_paired_tests.tsv
complexity/complexity_model.tsv
complexity/complexity_summary.txt
complexity/complexity_scores.tsv
complexity/*_stage_complexity.tsv
```

Environment variables:

| Variable | Default | Effect |
| --- | --- | --- |
| `CABLESIM_PCA_ORIENTATION` | `positive` | `positive` uses absolute PC1 loadings, `signed` keeps the signed projection |
| `CABLESIM_COMPLEXITY_FILTER` | `identified` | training rows: `identified`, `all`, `unidentified`, `non-unknot` |
| `CABLESIM_COMPLEXITY_PER_VARIANT` | `1` | `0` skips per-variant complexity output |

The driver also writes `catalog_with_complexity.tsv`,
`cut_eval_with_complexity.tsv`, and `cuts_to_zero_with_complexity.tsv` next to
the input tables.  These are consumed by `scripts/analyze_oracle.py`.

## Generate LaTeX tables

```bash
julia --project=cabletangleanalysis \
  cabletangleanalysis/scripts/make_latex_tables.jl \
  <dataset-root>/analysis_julia/<variant> <latex-out-dir>
```

Writes `table_pca.tex`, `table_complexity_first_cut.tex`,
`table_complexity_resolved.tex`, and `tables.tex`.

## Package API

Load the module from the project directory:

```julia
using Pkg; Pkg.activate("cabletangleanalysis")
using CableTangleAnalysis
```

All paths and options live in `Settings`; `DEFAULT_SETTINGS` is used when no
settings argument is passed.

```julia
settings = Settings(
    metrics_file = "path/to/catalog.tsv",
    cuts_file = "path/to/cuts_to_zero.tsv",
    cut_eval_file = "path/to/cut_eval.tsv",
    cut_eval_stage_file = "path/to/cut_eval_stage_metrics.tsv",
    cuts_to_zero_stage_file = "path/to/cuts_to_zero_stage_metrics.tsv",
)
```

Individual analyses:

```julia
plot_cuts_histogram(settings)              # cuts-per-strategy histogram
plot_cut_eval_comparison(settings)         # single-cut scatter vs heuristic
cuts_violin_comparisons(settings)          # four violins + paired tests

clear_data_of_unwanted_knots(settings)     # write cleaned PCA input
knot_distribution(settings)                # knot-type bar chart
run_pca(settings)                          # PCA figure and tables
correlation(settings)                      # correlation heatmap

model = analyze_complexity(settings)       # fit model, score stages, augment tables
plot_complexity_reduction_first_cut(settings; model = model)
plot_complexity_reduction_resolved(settings; model = model)
complexity_reduction_first_cut_paired_tests(settings; model = model)
complexity_reduction_resolved_paired_tests(settings; model = model)
```

`run_all(settings)` runs the same pipeline in one call.

## Project structure

```text
cabletangleanalysis/
|-- Project.toml
|-- README.md
|-- STATISTICS.md
|-- data/                  local input tables (not tracked)
|-- scripts/
|   |-- run_dataset.jl     dataset driver
|   |-- make_latex_tables.jl
|-- src/
|   |-- CableTangleAnalysis.jl  module, Settings, run_all
|   |-- data.jl                 loaders, knot grouping, cleaning
|   |-- utils.jl                TSV/column helpers and formatting
|   |-- plotutils.jl            label rendering and strategy colours
|   |-- paired.jl               paired tests, reports, violin plotting
|   |-- cuts.jl                 cut readers, histogram, cut comparisons, violins
|   |-- complexity.jl           PCA complexity model, stage scoring, violins
|   |-- pca.jl                  metric PCA and projection figure
|   |-- correlation.jl          correlation matrix, heatmap, redundancy check
|   |-- analysis.jl             knot-type distribution
```

## Data notes

- Undefined knot types (`none`, `unknown`, `??`, empty) are excluded from the
  metric, correlation, PCA, complexity, and cut analyses because their
  outside-knot geometry is undefined.  They are still counted in the knot
  distribution.
- Tangles that start with zero threading are excluded from the cut plots and
  paired tests since they cannot need a cut.
- Random strategies are averaged over their trials per seed; all strategy
  comparisons are paired by seed.
- The oracle is the shortest sequence within the hoop-midpoint action set, not
  a lower bound over arbitrary cable cuts.  See `STATISTICS.md`.
