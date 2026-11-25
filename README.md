# CableSim

CableSim is a Linux-focused fork of [PositionBasedDynamics](https://github.com/InteractiveComputerGraphics/PositionBasedDynamics) for interactive cable manipulation and knot/tanglement analysis. The main application, `MyCableDemo`, simulates one or more cable pieces with position-based rod constraints, self-collision, mouse dragging and cutting, and evaluates the resulting topology with a bundled copy of [KymoKnot](https://kymoknot.sissa.it/).

This branch (`linux`) contains the Linux port of the demo and the batch analysis workflows used to build and study the knot-seed catalog.

## Highlights

- Interactive cable simulation with grab-and-drag manipulation, cutting, end pulls, gravity, and ground contact.
- Random cable generation from reproducible seeds, plus mathematical knot presets such as `5_2`.
- Self-collision detection with a side-switch detector and rest-state creep.
- Tanglement analysis: hoop crossings, threading counts, knot spans, self-collision/close-pair statistics, curvature densities, and a catalog-level entanglement score.
- PCA-derived complexity score trained on non-threading geometric metrics.
- Knot type detection and localization through KymoKnot (top-down search by default).
- Unattended catalog generation, cut-heuristic evaluation, and cuts-to-zero evaluation.

## Requirements

- Linux (x86-64); the code uses POSIX paths such as `/proc/self/exe`.
- CMake 3.10 or newer.
- A C++17 compiler (the demo target explicitly requests C++17).
- [Eigen3](https://eigen.tuxfamily.org/).
- OpenGL and the usual X11/GLFW development packages.
- Git and network access during the first CMake configure, because CMake fetches pinned revisions of **Discregrid** and **GenericParameters**.
- Optional: OpenMP.
- Python 3 with NumPy and Matplotlib for complexity and oracle analysis (skip it with `--no-complexity` / `--no-oracle-analysis`).

Most other dependencies are bundled under `extern/` (GLFW, imgui, glew, freeglut, hapPLY, md5). KymoKnot is bundled under `KymoKnot-clone/`.

On Debian/Ubuntu, a typical dependency set is:

```bash
sudo apt install build-essential cmake git libeigen3-dev python3-numpy python3-matplotlib \
  libgl1-mesa-dev libx11-dev libxrandr-dev libxinerama-dev \
  libxcursor-dev libxi-dev
```

## Build

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"
```

The executable is written to `bin/MyCableDemo` (the project sets `CMAKE_RUNTIME_OUTPUT_DIRECTORY` to `bin/`).

If you do not need the Python bindings, disable them to keep the build smaller:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DUSE_PYTHON_BINDINGS=OFF
```

If Discregrid or GenericParameters are already installed, point CMake at them instead of fetching:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release \
  -DDiscregrid_INCLUDE_DIR=/path/to/include \
  -DDiscregrid_LIB=/path/to/libDiscregrid.a \
  -DDiscregrid_DEBUG_LIB=/path/to/libDiscregrid_d.a \
  -DGenericParameters_INCLUDE_DIR=/path/to/include
```

### KymoKnot

Knot-type analysis shells out to `KymoKnot-clone/KymoKnot_linear.x`. Build it once:

```bash
make -C KymoKnot-clone
```

If you use the binary already present in the repository, make sure it is executable:

```bash
chmod +x KymoKnot-clone/KymoKnot_linear.x
```

`MyCableDemo` looks for KymoKnot next to the executable (`bin/../KymoKnot-clone`) and falls back to `KymoKnot-clone` relative to the working directory. You can override the location:

```bash
export KYMOKNOT_ROOT="$PWD/KymoKnot-clone"
```

Without KymoKnot, the demo still runs, but knot-type reports are marked unavailable.

## Running

Start the interactive demo from the repository root:

```bash
./bin/MyCableDemo
```

Reproduce a specific random cable:

```bash
./bin/MyCableDemo --seed=42
```

Load a knot preset from a seed bank (the optional `:N` selects the Nth occurrence):

```bash
./bin/MyCableDemo --knot=5_2
./bin/MyCableDemo --knot=5_2:2
```

The seed bank is not checked in. It is created/appended in `Demos/MyCableDemo/` when a run detects a knot, or you can supply your own file with `--seed-bank=<path>`. A direct `--seed` takes precedence over `--knot`.

### Modes

| Mode | Flags | Output |
| --- | --- | --- |
| Interactive | *(default)* | Live simulation and analysis in the OpenGL window |
| Catalog generation | `--catalog-run` plus optional `--catalog-max`, `--catalog-start`, `--catalog-file` | Appends catalog rows to a TSV file |
| Cut heuristic evaluation | `--eval-run` plus optional `--eval-output`, `--eval-random-trials` | Replaces the cut-evaluation TSV |
| Cuts-to-zero evaluation | `--eval-cuts-to-zero` plus optional `--eval-ctz-output`, `--eval-ctz-max-cuts`, `--eval-ctz-random-trials` | Replaces the cuts-to-zero TSV |
| Topology audit | `--topology-audit-run` plus optional `--topology-audit-steps`, `--topology-audit-cuts`, `--topology-audit-max-seeds`, `--topology-audit-output` | Writes per-phase and per-event audit TSVs plus a summary |

The `--eval-*` and catalog modes run without user interaction, but they still create an OpenGL window. On a machine without a display, run them under a virtual display, for example:

```bash
xvfb-run -a ./bin/MyCableDemo --catalog-run --catalog-max=10
```

Catalog mode writes one sample after the rest-state creep has settled, then regenerates the cable for the next seed. When `--catalog-max` is reached, it logs the output path and exits automatically. Cut evaluation and cuts-to-zero read `knotType` and `seed` columns from the catalog file, so a catalog must exist first.

Example catalog run:

```bash
./bin/MyCableDemo \
  --catalog-run \
  --catalog-start=1000 \
  --catalog-max=50 \
  --catalog-file=Demos/MyCableDemo/knot_seed_catalog.tsv
```

Example cut-heuristic evaluation:

```bash
./bin/MyCableDemo \
  --catalog-file=Demos/MyCableDemo/knot_seed_catalog.tsv \
  --eval-run \
  --eval-output=Demos/MyCableDemo/knot_cut_eval.tsv \
  --eval-random-trials=10
```

Example cuts-to-zero evaluation:

```bash
./bin/MyCableDemo \
  --catalog-file=Demos/MyCableDemo/knot_seed_catalog.tsv \
  --eval-cuts-to-zero \
  --eval-ctz-output=Demos/MyCableDemo/knot_cuts_to_zero.tsv \
  --eval-ctz-max-cuts=10 \
  --eval-ctz-random-trials=20
```

### Topology audit

`--topology-audit-run` checks that the cable's topology is preserved while it settles under gravity and after each heuristic cut. It runs over the seeds in the catalog, picks the same hoop-midpoint cut used by the heuristic strategy, and audits **every simulation substep** after the initial settling phase and after every cut. There is deliberately no sampling stride: a sparse check could miss a pass-through.

A cut is an intentional topology change, so the auditor re-baselines immediately after the rebuilt geometry is in place. Only motion after the cut is judged.

```bash
xvfb-run -a ./bin/MyCableDemo \
  --topology-audit-run \
  --catalog-file=results/dataset_final/catalog.tsv \
  --topology-audit-max-seeds=100 \
  --topology-audit-steps=300 \
  --topology-audit-cuts=3 \
  --topology-audit-output=results/topology_audit
```

For a full catalog, [`scripts/run_topology_audit.sh`](scripts/run_topology_audit.sh) shards the catalog across parallel single-threaded workers and merges the per-shard results into one set of files:

```bash
scripts/run_topology_audit.sh \
  --catalog results/full_dataset_20261106_112203/current/catalog.tsv \
  --output-prefix results/full_dataset_20261106_112203/topology_audit \
  --workers 32 --steps 300 --cuts 3
```

Single-threaded workers are about three times faster per seed than letting OpenMP spread one seed over all cores, so sharding by seed is the efficient way to cover a large catalog. Per-shard catalogs, logs, and outputs are kept under `<prefix>_work/`.

| Argument | Description |
| --- | --- |
| `--topology-audit-run` | Enable unattended topology auditing over catalog seeds |
| `--topology-audit-steps=<n>` | Audited substeps per phase (default 1000) |
| `--topology-audit-cuts=<n>` | Heuristic cuts per seed (default 0) |
| `--topology-audit-max-seeds=<n>` | Limit the number of catalog seeds (default: all) |
| `--topology-audit-output=<prefix>` | Output prefix (default `<catalog-dir>/topology_audit`) |
| `--topology-audit-contact-eps=<value>` | Centre-line contact threshold (default `1e-4`) |
| `--topology-audit-trace` | Also write a per-step trace file |
| `--dt=<value>` | Physics timestep (default 0.005) |
| `--collision-iterations=<n>` | Self-collision solver iterations per substep (default 8) |
| `--collision-compliance=<value>` | Contact compliance (default 1e-10) |
| `--collision-max-correction=<value>` | Max per-endpoint contact correction (default 0.08) |
| `--normal-damping=<value>` | Contact normal damping (default 0.15) |
| `--tangential-damping=<value>` | Contact tangential damping, velocity-proportional (default 0.05); not Coulomb/stick–slip friction |
| `--hash-cell-size=<value>` | Self-collision spatial-hash cell size (default 0.2) |

#### How the audit works

Gravity here is along `-Y`, so the cable is a 2.5-D object: its over/under structure lives in the `XZ` projection, and `Y` records which strand is on top. That is why no full 3-D knot invariant is needed. For every pair of non-adjacent edges whose `XZ` projections cross, the auditor records:

```text
heightDelta = y_A - y_B   at the projected crossing point
```

Because both points share `x` and `z`, `|heightDelta|` is also the strands' centre-line distance at that crossing. For a ground-resting tangle, touching strands sit at `|heightDelta| ≈ 2 * ropeRadius`, so the default contact threshold is far below that and only real centre-line contact trips it.

Violation types:

| Type | Meaning |
| --- | --- |
| `CROSSING_FLIP` | A persistent crossing changed the sign of `heightDelta`: the strands passed through each other |
| `BIRTH_AT_CONTACT` | A crossing appeared while the strands were within `contactEpsilon`, so it was born through contact rather than by a projection rearrangement |
| `DEATH_AT_CONTACT` | A crossing disappeared while the strands were within `contactEpsilon` |
| `CENTERLINE_CONTACT` | A crossing stayed within `contactEpsilon`; the strands are touching or interpenetrating at the centre line |

Crossings that appear or disappear with real vertical separation are ordinary projection rearrangements (a Reidemeister-style move in the plan view) and are not reported. Adjacent edges of the same piece are excluded, and endpoint intersections are skipped to avoid double counting at shared particles.

#### Outputs

All files use the `--topology-audit-output` prefix:

| File | Contents |
| --- | --- |
| `<prefix>_seeds.tsv` | One row per seed and phase (`pre_cut`, `post_cut_1`, ...) |
| `<prefix>_events.tsv` | One row per violation |
| `<prefix>_summary.txt` | Aggregate totals/extrema and solver/config parameters across the run |
| `<prefix>_trace.tsv` | Optional per-step crossing, contact, and penetration diagnostics (`--topology-audit-trace`) |

`_seeds.tsv` records per-phase topology and solver-contact diagnostics, including:

```text
seed knotType phase steps cutApplied cutPiece cutEdge
violations centerlineContacts surfaceOverlaps crossingsPerStep minAbsHeight
candidatePairs penetratingContacts maxPenetration meanPenetration
residualCandidatePairs residualPenetratingContacts residualMaxPenetration residualMeanPenetration
maxGroundPenetration cappedCorrections correctionEndpoints maxAppliedCorrection normalFlips
firstViolationStep firstViolationType
firstPieceA firstEdgeA firstPieceB firstEdgeB
firstPreviousHeight firstCurrentHeight
```

`penetratingContacts`/`maxPenetration`/`meanPenetration` summarize pre-correction (detected) self-collision penetrations accumulated over solver iterations in audited substeps. `residual*` fields are measured after the final self-collision iteration of each substep, before the ground clamp, and `maxGroundPenetration` records pre-clamp ground penetration. These diagnostics separate *detected-and-corrected* interpenetration from *remaining* penetration/contact-violation residuals at audited states.

`_events.tsv` includes the crossing/topology event plus the local contact-penetration diagnostics at the event step:

```text
seed knotType phase step type
pieceA edgeA pieceB edgeB x z
previousX previousZ previousHeight currentHeight
candidatePairs penetratingContacts maxPenetration
residualPenetratingContacts residualMaxPenetration
maxGroundPenetration normalFlips
```

`surfaceOverlaps` counts crossings with `|heightDelta| < 2 * ropeRadius`. It is a geometric diagnostic only: in a ground-resting tangle most crossings are at that distance, so it does not by itself indicate a topology violation.

Read the results as **observed discrete-substep evidence, not a formal proof**. The simulation advances with finite `dt`, a finite number of compliance-softened, capped/averaged collision-correction passes per substep, and finite contact passes. Every-substep auditing observes every configured substep after settling/re-baselining, but a crossing or penetration event that is born and resolved entirely between audited states would not be seen. The honest claim is: *no topology violation and no residual penetration above the reported audited-state diagnostics were observed at audited substeps under the stated solver configuration.* For the mathematical formulation and sensitivity methodology, see [`docs/topology_audit_math.md`](docs/topology_audit_math.md) and [`docs/creep_sensitivity_math.md`](docs/creep_sensitivity_math.md).

### Full parallel pipeline

[`scripts/run_full_pipeline.sh`](scripts/run_full_pipeline.sh) runs the catalog and both evaluations for a seed range, split across parallel workers:

```bash
scripts/run_full_pipeline.sh --start 1 --end 1000 --workers 8
```

By default it generates three rest-state-creep variants:

| Directory | Creep value |
| --- | ---: |
| `zero/` | 0 |
| `current/` | 0.02 |
| `high/` | 0.2 |

Each variant contains `catalog.tsv`, `cut_eval.tsv`, `cuts_to_zero.tsv`, `cut_eval_stage_metrics.tsv`, `cuts_to_zero_stage_metrics.tsv`, a `summary.txt`, and per-worker logs. The stage-metrics files contain the five complexity metrics for every tangle before and after each cut, so the PCA score can be recomputed for every stage independently. A top-level `manifest.tsv` lists all variants. Use `--creep <value>` to run a single creep value instead (written directly to `--out`), or override the high/current values with `--creep-high` / `--creep-current`. Per-variant complexity analysis runs for every variant; the pooled model and oracle analysis run in the three-variant mode.

The script also sets `--creep-settle-steps=300` and lowers the creep activation threshold to `0.0001`. With the code default threshold (`0.01`), sustained mismatch under gravity-only settling stays below the threshold, so the creep variants would otherwise be identical. Override these with `--creep-settle-steps` and `--creep-threshold`.

The script sets `OMP_NUM_THREADS=1` for each worker and starts a virtual X display, because batch mode still creates an OpenGL window. It runs the exhaustive single-cut oracle by default; pass `--no-oracle` to skip it. Note that seed `0` is normalized to `1` by `MyCableDemo`; use `--start 1 --end 1000` for 1000 unique seeds. Running all three creep variants takes roughly three times as long as a single variant.

After the simulation stages it runs the Julia analysis in `cabletangleanalysis/`, which produces the knot-type distribution, PCA plot, correlation heatmap, cut-evaluation figures, and the complexity analysis under `analysis_julia/<variant>/`. The pooled complexity model is written to `complexity/`, and each creep variant gets its own model under `analysis_julia/<variant>/complexity/`. Set `JULIA_BIN` or `JULIA_PROJECT` if Julia is not on `PATH` or the project is elsewhere. Pass `--no-complexity` to skip the Julia pass, or `--no-complexity-per-variant` to fit only the pooled model.

Finally, [`scripts/analyze_oracle.py`](scripts/analyze_oracle.py) analyses oracle single-cut performance against threading count and complexity score under `analysis_oracle/`. It requires Matplotlib; skip it with `--no-oracle-analysis`.

### Cut strategies

The evaluation modes compare these cut strategies:

| Strategy | Selection | Modes |
| --- | --- | --- |
| Heuristic | Hoop with the highest `pairSumScore` | both |
| Random edge | Uniformly random valid cable edge | both |
| Random hoop | Uniformly random hoop with a positive score | both |
| Biggest hoop | Hoop with the largest edge span (`endEdge - startEdge`) | both |
| Smallest hoop | Hoop with the smallest edge span | both |
| Max curvature | Valid edge with the largest discrete bend angle | both |
| Oracle | Every valid edge; best single-cut reduction | cut evaluation, opt-in |

Hoop-based strategies cut at the hoop midpoint. Max-curvature cuts use the largest turning angle at the edge endpoints and respect `excludeArcLength`. Random strategies use a per-seed RNG, so their results are reproducible. `--eval-random-trials` controls both random-edge and random-hoop trials in cut evaluation; `--eval-ctz-random-trials` controls both for cuts-to-zero.

`MyCableDemo` computes the exhaustive oracle only when `--eval-oracle` is passed. The pipeline script enables it by default; disable it with `--no-oracle`.

Cut evaluation appends these columns:
`delta_random_hoop_mean`, `delta_random_hoop_std`, `delta_biggest_hoop`, `biggest_hoop_span`, `biggest_hoop_score`, `hoop_count`, `delta_smallest_hoop`, `smallest_hoop_span`, `smallest_hoop_score`, `delta_max_curvature`, `max_curvature`, `delta_oracle`, `oracle_evaluated`, `oracle_piece_index`, `oracle_edge_index`.

Cuts-to-zero appends:
`random_hoop_cuts_mean`, `random_hoop_cuts_std`, `random_hoop_resolved_rate`, `biggest_hoop_cuts`, `biggest_hoop_resolved`, `hoop_count`, `smallest_hoop_cuts`, `smallest_hoop_resolved`, `max_curvature_cuts`, `max_curvature_resolved`.

### Complexity score

The complexity score is computed after catalog generation by the Julia module in `cabletangleanalysis/`. It fits a PCA on five geometric metrics that do not directly determine the threading count:

```text
segmentCount
selfCollisionCount
closePairCount
totalCurvature
outsideKnotCurvature
```

The first principal component is used as the metric weighting. By default (`--pca-orientation positive`), the absolute PC1 loadings are used so every metric contributes in the positive complexity direction: higher values always mean a more complex tangle. The model file writes the actual weight in `pc1_loading` and the raw signed PC1 loading in `pc1_loading_signed`.

Use `--pca-orientation signed` to keep the exact signed PC1 projection instead. In that mode individual loadings may be negative because PC1 describes a trade-off direction.

By default, PCA is fitted on rows identified by KymoKnot (`--complexity-filter identified`), excluding `none` and `??` rows. If fewer than three rows match, the analysis falls back to all rows and prints a warning. Other filters are `all`, `unidentified`, and `non-unknot`. When the pipeline runs multiple creep variants, a pooled model is fitted across all variants, and a separate model is fitted per creep variant so each setting has its own weights and scores.

`MyCableDemo` writes `cut_eval_stage_metrics.tsv` and `cuts_to_zero_stage_metrics.tsv` with the five complexity metrics for the initial tangle and for every cut stage. Julia applies the fitted model to those rows, so complexity reduction can be compared across strategies without recomputing the PCA during simulation.

Outputs:

| File | Contents |
| --- | --- |
| `complexity/complexity_model.tsv` | Pooled feature means, standard deviations, oriented PC1 weight, and signed PC1 loading |
| `complexity/complexity_summary.txt` | Explained variance and orientation mode |
| `complexity/complexity_scores.tsv` | Per-seed weighted score, signed PC1, normalized score, percentile, raw features |
| `analysis_julia/<variant>/complexity/complexity_model.tsv` | Per-variant PCA model |
| `analysis_julia/<variant>/complexity/cut_eval_stage_complexity.tsv` | Single-cut stage metrics with `delta_threading` and `delta_complexity` |
| `analysis_julia/<variant>/complexity/cuts_to_zero_stage_complexity.tsv` | All cuts-to-zero stages with both deltas |
| `<variant>/catalog_with_complexity.tsv` | Catalog with `complexity_pc1`, `complexity_pc1_signed`, `complexity_score`, `complexity_percentile` |
| `<variant>/cut_eval_with_complexity.tsv` | Cut evaluation joined to the initial complexity score |
| `<variant>/cuts_to_zero_with_complexity.tsv` | Cuts-to-zero evaluation joined to the initial complexity score |
| `analysis_julia/<variant>/figures/cut_eval/cut_eval_comparison.pdf` | Single-cut comparison using threading reduction |
| `analysis_julia/<variant>/figures/cut_eval/cut_eval_comparison_complexity.pdf` | Single-cut comparison using complexity reduction |

The simulation and strategy evaluation remain complexity-agnostic; the PCA score is applied in post-processing. The legacy `entanglementScore` column is still written by the catalog, but `complexity_score` is the PCA-derived score.

### Oracle analysis

`scripts/analyze_oracle.py` is the second-pass oracle analysis. It uses the per-variant complexity scores from the first pass and the `delta_oracle` results from cut evaluation, and reports:

- oracle threading removal by initial threading count `T_before`
- oracle threading removal by complexity quartile
- oracle fraction removed (`delta_oracle / T_before`) by both dimensions
- how often the oracle beats the heuristic strategy

Outputs per variant under `analysis_oracle/<variant>/`:

```text
oracle_summary.tsv
oracle_by_threading.tsv
oracle_by_complexity.tsv
oracle_by_threading_complexity.tsv
oracle_by_threading.pdf
oracle_by_complexity.pdf
oracle_heatmap.pdf
```

To recompute the complexity and cut-evaluation analysis on an existing dataset:

```bash
julia --project=cabletangleanalysis \
  cabletangleanalysis/scripts/run_dataset.jl \
  results/full_dataset_<timestamp> zero current high
```

or from the pipeline:

```bash
CABLESIM_PCA_ORIENTATION=positive \
CABLESIM_COMPLEXITY_FILTER=identified \
julia --project=cabletangleanalysis \
  cabletangleanalysis/scripts/run_dataset.jl results/full_dataset_<timestamp>
```

For the reduced analysis set (knot distribution, PCA, correlation, cut-eval
comparisons, cuts histogram, and PCA tables), run the same script; the
figure set is fixed and no separate core driver exists any more.

The outputs live under `analysis_julia/<variant>/` and include
`tables/pca_results.tsv`, `tables/pca_scores.tsv`,
`figures/cut_eval/cut_eval_comparison_complexity.pdf` (first-cut complexity
reduction), and `figures/cut_eval/complexity_reduction_resolved.pdf`
(complexity reduction after threading resolution). Both complexity plots show
the signed PC1 reduction.

### Command-line reference

| Argument | Description |
| --- | --- |
| `--no-cache` | Disable the framework's internal cache |
| `--seed-bank=<path>` | Knot-seed preset file (default: `Demos/MyCableDemo/knot_seed_bank.txt`) |
| `--catalog-file=<path>` | Catalog TSV for catalog output and evaluation input (default: `Demos/MyCableDemo/knot_seed_catalog.tsv`) |
| `--seed=<n>`, `--random-seed=<n>` | Fixed RNG seed for the random cable |
| `--knot=<name>`, `--knot-type=<name>`, `--knot <name>` | Load a named seed preset; `name:N` selects the Nth occurrence |
| `--catalog-run` | Enable catalog mode |
| `--catalog-render` | Keep rendering in catalog mode |
| `--catalog-steps-per-frame=<n>` | Simulation steps per frame in catalog mode; also enables catalog mode |
| `--catalog-max=<n>` | Stop after N catalog samples; also enables catalog mode |
| `--catalog-start=<n>` | Start the catalog seed sequence at N; also enables catalog mode |
| `--creep=<value>` | Rest-state creep rate $\gamma$; `0` disables adaptation but keeps batch triggering enabled |
| `--creep-settle-steps=<n>` | Steps to run after creep activation/capture before batch evaluation (default 0) |
| `--creep-threshold=<value>` | Creep mismatch activation threshold (default 0.01; pipelines/audits commonly use 0.0001 under gravity-only settling) |
| `--eval-run` | Run cut evaluation (all non-oracle strategies) |
| `--eval-oracle` | Also run the exhaustive single-cut oracle; implies `--eval-run` |
| `--eval-cuts-to-zero` | Run cuts-to-zero evaluation |
| `--eval-output=<path>` | Output file for cut-heuristic evaluation |
| `--eval-stage-metrics-output=<path>` | Per-stage complexity metrics for cut evaluation (default: derived from `--eval-output`) |
| `--eval-ctz-output=<path>` | Output file for cuts-to-zero evaluation |
| `--eval-ctz-stage-metrics-output=<path>` | Per-stage complexity metrics for cuts-to-zero (default: derived from `--eval-ctz-output`) |
| `--eval-random-trials=<n>` | Random edge/hoop trials per sample in cut evaluation (default 10) |
| `--eval-ctz-max-cuts=<n>` | Maximum cuts per sample in cuts-to-zero evaluation (default 10) |
| `--eval-ctz-random-trials=<n>` | Random edge/hoop trials per sample in cuts-to-zero evaluation (default 20) |
| `--dt=<value>` | Physics timestep (default 0.005) |
| `--collision-iterations=<n>` | Self-collision solver iterations per substep (default 8) |
| `--collision-compliance=<value>` | Contact compliance (default 1e-10) |
| `--collision-max-correction=<value>` | Max per-endpoint contact correction (default 0.08) |
| `--normal-damping=<value>` | Contact normal damping (default 0.15) |
| `--tangential-damping=<value>` | Contact tangential damping, velocity-proportional (default 0.05); not Coulomb/stick–slip friction |
| `--hash-cell-size=<value>` | Self-collision spatial-hash cell size (default 0.2) |
| `--topology-audit-run` | Run unattended topology auditing over catalog seeds |
| `--topology-audit-steps=<n>` | Audited substeps per phase (default 1000) |
| `--topology-audit-cuts=<n>` | Heuristic cuts per seed (default 0) |
| `--topology-audit-max-seeds=<n>` | Limit the number of catalog seeds (default: all) |
| `--topology-audit-output=<prefix>` | Output prefix (default `<catalog-dir>/topology_audit`) |
| `--topology-audit-contact-eps=<value>` | Centre-line contact threshold (default `1e-4`) |
| `--topology-audit-trace` | Also write a per-step trace file |

Default output paths are resolved relative to the executable, so running `bin/MyCableDemo` uses the checked-in `Demos/MyCableDemo/` directory. The `*_v3.tsv` files there are reference artifacts from earlier runs; the defaults use the unsuffixed names unless you pass an explicit path.

## Controls

### Mouse

| Input | Action |
| --- | --- |
| Left-click + drag on a cable | Grab and drag the selected cable segment |
| Left-click on a hoop marker | Select the hoop and print its tangle score |
| Left-drag on empty space with `Alt` | Rotate the view |
| Left-drag on empty space with `Shift` | Pan the view |
| Left-drag on empty space with `Ctrl` | Dolly the view |
| Mouse wheel | Zoom |

Clicks on cables and hoop markers are handled by the demo; clicks on empty space fall through to the camera controls.

### Keyboard

| Key | Action |
| --- | --- |
| `r` | Reset the current simulation |
| `o` | Regenerate a random cable with a new seed |
| `k` | Toggle the mathematical knot preset (`5_2` first, then a straight cable) |
| `t` | Run and print tanglement analysis for the current cable |
| `c` | Cut the last selected cable segment |
| `v` | While dragging, enable upward pull on the grabbed segment |
| `g` | Toggle gravity |
| `*` | Toggle end pull forces |
| `+` / `-` | Increase / decrease end pull strength |
| `z` / `x` | Move the camera toward / away from the scene |
| Arrow keys | Pan the camera |
| `1`–`4` | Rotate the view around the x / y axes |
| `Esc` | Quit |

## Analysis and data

Pressing `t` prints a report containing the detected knot type(s), hoop/threading metrics, self-collision statistics, knot span, and proximity/curvature scores. Hoop markers are generated from the analysis and can be clicked for per-hoop scores.

Catalog rows include, among others:

- `knotType`, `seed`, `pieceIndex`, `knotCount`, `pieceCount`
- hoop and threading counts: `totalHoopCount`, `threadedHoopCount`, `totalThreadingCount`
- collision statistics: `selfCollisionCount`, `closePairCount`, `contactParticipationCount`, `minPairDistance`
- topology/geometry: `knotSpanFraction`, `proximityScore`, `totalCurvature`, `outsideKnotCurvature`, `totalCurvatureDensity`, `outsideKnotCurvatureDensity`
- `entanglementScore`, a catalog-normalized combination of threading, close pairs, contact participation, knot span, outside-knot curvature, extra knots, and proximity

## Configuration

Most simulation defaults live in [`Demos/MyCableDemo/CablePhysics.h`](Demos/MyCableDemo/CablePhysics.h):

- `RopeConfig`: rod stiffness, segment length, render radius
- `PhysicsConfig`: time step, end-pull spring/damping, velocity blending
- `RandomCableConfig`: control-point count, lateral bounds, smoothing, seed
- `SelfCollisionConfig`: radius, compliance, iterations, damping, hash cell size
- `RestShapeCreepConfig`: creep rate, activation delay, evaluation delay, update interval, and activation threshold

Scene construction, keyboard bindings, evaluation orchestration, and output paths are in [`Demos/MyCableDemo/main.cpp`](Demos/MyCableDemo/main.cpp).

## Repository layout

| Path | Contents |
| --- | --- |
| `Demos/MyCableDemo/` | Cable demo, physics helpers, analyzers |
| `scripts/` | Dataset pipeline, parallel topology audit runner, and PCA complexity analysis |
| `Demos/Visualization/` | OpenGL/imgui visualization helpers |
| `PositionBasedDynamics/` | PBD/XPBD constraint library |
| `Simulation/` | Simulation framework |
| `Utils/` | File system, scene loading, and utility code |
| `extern/` | Bundled third-party libraries |
| `KymoKnot-clone/` | Bundled KymoKnot knot analyzer |
| `cabletangleanalysis/` | Julia analysis submodule (PCA, correlations, cut-strategy comparisons) |
| `data/` | Runtime shaders (geometry is generated procedurally) |
| `bin/` | Build output (ignored by Git) |
| `build/` | CMake build tree (ignored by Git) |

The Julia analysis lives in the `cabletangleanalysis` submodule. After cloning, initialize it with:

```bash
git submodule update --init --recursive
```

## Troubleshooting

- **KymoKnot reports unavailable**: build `KymoKnot-clone` with `make -C KymoKnot-clone`, ensure `KymoKnot_linear.x` is executable, or set `KYMOKNOT_ROOT`.
- **CMake configure fails while cloning Discregrid/GenericParameters**: check Git/network access or provide the preinstalled paths shown above.
- **Evaluation mode fails to open a window**: unattended modes still need an OpenGL context; use `xvfb-run -a`.
- **Seed bank or catalog not found**: pass `--seed-bank` / `--catalog-file`, or run from the repository root so the default paths resolve.
- **Complexity analysis fails with Julia missing**: install Julia and instantiate the `cabletangleanalysis` project (`julia --project=cabletangleanalysis -e 'using Pkg; Pkg.instantiate()'`), set `JULIA_BIN` / `JULIA_PROJECT`, or pass `--no-complexity`.
- **Oracle analysis fails with missing Python packages**: install `requirements.txt` in a virtual environment, or pass `--no-oracle-analysis`.
- **Knot type looks wrong for complex knots**: KymoKnot uses Alexander determinants at `t=-1` and `t=-2`; complex knots can share determinants with other knots. See `KymoKnot-clone/README.md` for details.

## License and upstream

CableSim is distributed under the MIT license; see [`LICENSE`](LICENSE). It is based on PositionBasedDynamics by Jan Bender and contributors. KymoKnot is a separate package with its own license in [`KymoKnot-clone/COPYING.txt`](KymoKnot-clone/COPYING.txt).

Selected references:

- J. Bender, M. Müller, M. Macklin, "Position-Based Simulation Methods in Computer Graphics", Eurographics Tutorial, 2015.
- C. Deul, T. Kugelstadt, M. Weiler, J. Bender, "Direct Position-Based Solver for Stiff Rods", Computer Graphics Forum, 2018.
- T. Kugelstadt, E. Schoemer, "Position and Orientation Based Cosserat Rods", SCA, 2016.
- L. Tubiana, E. Orlandini, C. Micheletti, "Probing the Entanglement and Locating Knots in Ring Polymers", Progress of Theoretical Physics Supplement, 2011.
