# CableSim

CableSim is a Linux-focused fork of [PositionBasedDynamics](https://github.com/InteractiveComputerGraphics/PositionBasedDynamics) for interactive cable manipulation and knot/tanglement analysis. The main application, `MyCableDemo`, simulates one or more cable pieces with position-based rod constraints, self-collision, mouse dragging and cutting, and evaluates the resulting topology with a bundled copy of [KymoKnot](https://kymoknot.sissa.it/).

This branch (`linux`) contains the Linux port of the demo and the batch analysis workflows used to build and study the knot-seed catalog.

## Highlights

- Interactive cable simulation with grab-and-drag manipulation, cutting, end pulls, gravity, and ground contact.
- Random cable generation from reproducible seeds, plus mathematical knot presets such as `5_2`.
- Self-collision detection with a side-switch detector and rest-state creep.
- Tanglement analysis: hoop crossings, threading counts, knot spans, self-collision/close-pair statistics, curvature densities, and a catalog-level entanglement score.
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
- Optional: OpenMP, Python 3 with NumPy for the PCA script.

Most other dependencies are bundled under `extern/` (GLFW, imgui, glew, freeglut, hapPLY, md5). KymoKnot is bundled under `KymoKnot-clone/`.

On Debian/Ubuntu, a typical dependency set is:

```bash
sudo apt install build-essential cmake git libeigen3-dev \
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

Results are written to `results/catalog_<timestamp>/` as `catalog.tsv`, `cut_eval.tsv`, and `cuts_to_zero.tsv`, with per-worker logs under `logs/`. The script sets `OMP_NUM_THREADS=1` for each worker and starts a virtual X display, because batch mode still creates an OpenGL window.

Use `--help` for all options. Note that seed `0` is normalized to `1` by `MyCableDemo`; use `--start 1 --end 1000` for 1000 unique seeds.

### Cut strategies

The evaluation modes compare four cut strategies:

| Strategy | Selection |
| --- | --- |
| Heuristic | Hoop with the highest `pairSumScore` |
| Random edge | Uniformly random valid cable edge (existing baseline) |
| Random hoop | Uniformly random hoop with a positive score |
| Biggest hoop | Hoop with the largest edge span (`endEdge - startEdge`), tie-broken by score |

Hoop-based strategies cut at the hoop midpoint. Random strategies use a per-seed RNG, so their results are reproducible. `--eval-random-trials` controls both random-edge and random-hoop trials in cut evaluation; `--eval-ctz-random-trials` controls both for cuts-to-zero.

Cut evaluation appends these columns:
`delta_random_hoop_mean`, `delta_random_hoop_std`, `delta_biggest_hoop`, `biggest_hoop_span`, `biggest_hoop_score`, `hoop_count`.

Cuts-to-zero appends:
`random_hoop_cuts_mean`, `random_hoop_cuts_std`, `random_hoop_resolved_rate`, `biggest_hoop_cuts`, `biggest_hoop_resolved`, `hoop_count`.

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
| `--eval-run` | Run cut-heuristic evaluation |
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
- topology/geometry: `knotSpanFraction`, `proximityScore`, `totalCurvatureDensity`, `outsideKnotCurvatureDensity`
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
- `RestShapeCreepConfig`: rest-shape creep activation and rate

Scene construction, keyboard bindings, evaluation orchestration, and output paths are in [`Demos/MyCableDemo/main.cpp`](Demos/MyCableDemo/main.cpp).

## Repository layout

| Path | Contents |
| --- | --- |
| `Demos/MyCableDemo/` | Cable demo, physics helpers, analyzers, evaluation data, PCA script |
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
- **Knot type looks wrong for complex knots**: KymoKnot uses Alexander determinants at `t=-1` and `t=-2`; complex knots can share determinants with other knots. See `KymoKnot-clone/README.md` for details.

## License and upstream

CableSim is distributed under the MIT license; see [`LICENSE`](LICENSE). It is based on PositionBasedDynamics by Jan Bender and contributors. KymoKnot is a separate package with its own license in [`KymoKnot-clone/COPYING.txt`](KymoKnot-clone/COPYING.txt).

Selected references:

- J. Bender, M. Müller, M. Macklin, "Position-Based Simulation Methods in Computer Graphics", Eurographics Tutorial, 2015.
- C. Deul, T. Kugelstadt, M. Weiler, J. Bender, "Direct Position-Based Solver for Stiff Rods", Computer Graphics Forum, 2018.
- T. Kugelstadt, E. Schoemer, "Position and Orientation Based Cosserat Rods", SCA, 2016.
- L. Tubiana, E. Orlandini, C. Micheletti, "Probing the Entanglement and Locating Knots in Ring Polymers", Progress of Theoretical Physics Supplement, 2011.
