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
- PCA helper script for exploring catalog metrics and generating an HTML scatter plot.

## Requirements

- Linux (x86-64); the code uses POSIX paths such as `/proc/self/exe`.
- CMake 3.10 or newer.
- A C++17 compiler (the demo target explicitly requests C++17).
- [Eigen3](https://eigen.tuxfamily.org/).
- OpenGL and the usual X11/GLFW development packages.
- Git and network access during the first CMake configure, because CMake fetches pinned revisions of **Discregrid** and **GenericParameters**.
- Optional: OpenMP.
- Python 3 with NumPy for the PCA complexity pass (skip it with `--no-complexity`).

Most other dependencies are bundled under `extern/` (GLFW, imgui, glew, freeglut, hapPLY, md5). KymoKnot is bundled under `KymoKnot-clone/`.

On Debian/Ubuntu, a typical dependency set is:

```bash
sudo apt install build-essential cmake git libeigen3-dev python3-numpy \
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

Each variant contains `catalog.tsv`, `cut_eval.tsv`, `cuts_to_zero.tsv`, a `summary.txt`, and per-worker logs. A top-level `manifest.tsv` lists all variants. Use `--creep <value>` to run a single creep value instead (written directly to `--out`), or override the high/current values with `--creep-high` / `--creep-current`.

The script also sets `--creep-settle-steps=300` and lowers the creep activation threshold to `0.0001`. With the code default threshold (`0.01`), sustained mismatch under gravity-only settling stays below the threshold, so the creep variants would otherwise be identical. Override these with `--creep-settle-steps` and `--creep-threshold`.

The script sets `OMP_NUM_THREADS=1` for each worker and starts a virtual X display, because batch mode still creates an OpenGL window. It runs the exhaustive single-cut oracle by default; pass `--no-oracle` to skip it. Note that seed `0` is normalized to `1` by `MyCableDemo`; use `--start 1 --end 1000` for 1000 unique seeds. Running all three creep variants takes roughly three times as long as a single variant.

After the simulation stages it runs [`scripts/compute_complexity.py`](scripts/compute_complexity.py), which fits the PCA complexity model and joins the resulting score into the catalog and strategy outputs under `complexity/`. This requires Python 3 with NumPy. Pass `--no-complexity` to skip it.

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

The complexity score is computed after catalog generation by `scripts/compute_complexity.py`. It fits a PCA on five geometric metrics that do not directly determine the threading count:

```text
segmentCount
selfCollisionCount
closePairCount
totalCurvature
outsideKnotCurvature
```

The first principal component is used as the metric weighting. The score is the PC1 projection of the standardized metrics; the loadings are written to `complexity_model.tsv`. PC1 is sign-oriented so the summed loadings are positive, but individual loadings may be negative because PC1 describes a trade-off direction.

By default, PCA is fitted on rows identified by KymoKnot (`--complexity-filter identified`), excluding `none` and `??` rows. If fewer than three rows match, the script falls back to all rows and prints a warning. Other filters are `all`, `unidentified`, and `non-unknot`. When the pipeline runs multiple creep variants, all variants are pooled so the complexity score uses one common scale. To fit a separate model per variant, call `compute_complexity.py` once per catalog instead.

Outputs:

| File | Contents |
| --- | --- |
| `complexity/complexity_model.tsv` | Feature means, standard deviations, and PC1 weights |
| `complexity/complexity_summary.txt` | Training rows, filter, explained variance, orientation |
| `complexity/complexity_scores.tsv` | Per-seed PC1, normalized score, percentile, raw features |
| `<variant>/catalog_with_complexity.tsv` | Catalog with `complexity_pc1`, `complexity_score`, `complexity_percentile` |
| `<variant>/cut_eval_with_complexity.tsv` | Cut evaluation joined to the complexity score |
| `<variant>/cuts_to_zero_with_complexity.tsv` | Cuts-to-zero evaluation joined to the complexity score |
| `complexity/strategy_by_complexity.tsv` | Strategy means/medians by complexity quartile |
| `complexity/strategy_correlations.tsv` | Pearson correlation between complexity and strategy metrics |

The simulation and strategy evaluation remain complexity-agnostic; the PCA score is applied in post-processing. The legacy `entanglementScore` column is still written by the catalog, but `complexity_score` is the PCA-derived score.

To recompute the complexity pass on an existing dataset:

```bash
python3 scripts/compute_complexity.py \
  results/full_dataset_<timestamp>/{zero,current,high}/catalog.tsv \
  --out-dir results/full_dataset_<timestamp>/complexity
```

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
| `--creep=<value>` | Rest-state creep rate; `0` disables adaptation but keeps batch triggering enabled |
| `--creep-settle-steps=<n>` | Steps to run after creep activation before batch evaluation (default 0) |
| `--creep-threshold=<value>` | Creep mismatch activation threshold (default 0.01) |
| `--eval-run` | Run cut evaluation (all non-oracle strategies) |
| `--eval-oracle` | Also run the exhaustive single-cut oracle; implies `--eval-run` |
| `--eval-cuts-to-zero` | Run cuts-to-zero evaluation |
| `--eval-output=<path>` | Output file for cut-heuristic evaluation |
| `--eval-ctz-output=<path>` | Output file for cuts-to-zero evaluation |
| `--eval-random-trials=<n>` | Random edge/hoop trials per sample in cut evaluation (default 10) |
| `--eval-ctz-max-cuts=<n>` | Maximum cuts per sample in cuts-to-zero evaluation (default 10) |
| `--eval-ctz-random-trials=<n>` | Random edge/hoop trials per sample in cuts-to-zero evaluation (default 20) |

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

Reference artifacts in `Demos/MyCableDemo/`:

| File | Purpose |
| --- | --- |
| `knot_seed_catalog_v3.tsv` | Full catalog example |
| `knot_seed_catalog_v3_sorted_by_type_crossings.tsv` | Catalog sorted by knot type and crossing count |
| `knot_seed_catalog_v3_trimmed_metrics.tsv` | Reduced metric columns |
| `knot_seed_catalog_v3_pca_A_density.tsv` | Density-normalized PCA input |
| `knot_seed_catalog_v3_pca_A_density_pca_*` | PCA coordinates, loadings, and plot |
| `knot_cut_eval_v3.tsv` | Cut-heuristic evaluation example |
| `knot_cuts_to_zero_v3.tsv` | Cuts-to-zero evaluation example |

The checked-in `knot_cut_eval_v3.tsv` and `knot_cuts_to_zero_v3.tsv` files predate the random-hoop and biggest-hoop strategies; new evaluation runs include the additional columns.

### PCA workflow

`run_pca.py` is the older PCA visualization helper and defaults to density/threading features. The dataset pipeline's complexity score uses `scripts/compute_complexity.py` with the five non-threading metrics described above.

After producing a density-normalized metrics table, run:

```bash
python3 Demos/MyCableDemo/run_pca.py \
  Demos/MyCableDemo/knot_seed_catalog_v3_pca_A_density.tsv
```

The script requires NumPy and writes `<prefix>_pca_coords.tsv`, `<prefix>_pca_loadings.tsv`, and `<prefix>_pca_plot.html`. Useful options:

```bash
python3 Demos/MyCableDemo/run_pca.py INPUT.tsv \
  --features totalHoopDensity threadedHoopDensity threadingDensity \
  --pc-count 3 \
  --output-prefix my_analysis \
  --title "My cable PCA" \
  --spread 0.25
```

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
| `Demos/MyCableDemo/` | Cable demo, physics helpers, analyzers, evaluation data, legacy PCA script |
| `scripts/` | Dataset pipeline and PCA complexity analysis |
| `Demos/Visualization/` | OpenGL/imgui visualization helpers |
| `PositionBasedDynamics/` | PBD/XPBD constraint library |
| `Simulation/` | Simulation framework |
| `Utils/` | File system, scene loading, and utility code |
| `extern/` | Bundled third-party libraries |
| `KymoKnot-clone/` | Bundled KymoKnot knot analyzer |
| `data/` | Shaders, models, and scene data |
| `bin/` | Build output (ignored by Git) |
| `build/` | CMake build tree (ignored by Git) |

## Troubleshooting

- **KymoKnot reports unavailable**: build `KymoKnot-clone` with `make -C KymoKnot-clone`, ensure `KymoKnot_linear.x` is executable, or set `KYMOKNOT_ROOT`.
- **CMake configure fails while cloning Discregrid/GenericParameters**: check Git/network access or provide the preinstalled paths shown above.
- **Evaluation mode fails to open a window**: unattended modes still need an OpenGL context; use `xvfb-run -a`.
- **Seed bank or catalog not found**: pass `--seed-bank` / `--catalog-file`, or run from the repository root so the default paths resolve.
- **Complexity pass fails with NumPy missing**: install `python3-numpy` or pass `--no-complexity` to the pipeline.
- **Knot type looks wrong for complex knots**: KymoKnot uses Alexander determinants at `t=-1` and `t=-2`; complex knots can share determinants with other knots. See `KymoKnot-clone/README.md` for details.

## License and upstream

CableSim is distributed under the MIT license; see [`LICENSE`](LICENSE). It is based on PositionBasedDynamics by Jan Bender and contributors. KymoKnot is a separate package with its own license in [`KymoKnot-clone/COPYING.txt`](KymoKnot-clone/COPYING.txt).

Selected references:

- J. Bender, M. Müller, M. Macklin, "Position-Based Simulation Methods in Computer Graphics", Eurographics Tutorial, 2015.
- C. Deul, T. Kugelstadt, M. Weiler, J. Bender, "Direct Position-Based Solver for Stiff Rods", Computer Graphics Forum, 2018.
- T. Kugelstadt, E. Schoemer, "Position and Orientation Based Cosserat Rods", SCA, 2016.
- L. Tubiana, E. Orlandini, C. Micheletti, "Probing the Entanglement and Locating Knots in Ring Polymers", Progress of Theoretical Physics Supplement, 2011.
