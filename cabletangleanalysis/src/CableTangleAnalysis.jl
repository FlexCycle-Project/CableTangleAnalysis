module CableTangleAnalysis

# ------------------------------------------------------------------
# CableTangleAnalysis.jl — Main module entry point
#
# Analysis of 3D-printed cable tangle metrics: PCA, correlation,
# cut-evaluation, knot distribution, and a PCA-derived complexity score.
#
# Usage (from the project directory):
#   julia> using Pkg; Pkg.activate(".")
#   julia> using CableTangleAnalysis
#   julia> CableTangleAnalysis.run_all()
#
# The dataset driver used for the paper is scripts/run_dataset.jl.
# ------------------------------------------------------------------

# ---- Configuration struct ------------------------------------------------

"""
    Settings

Central configuration for the analysis pipeline.  Create a custom instance
to change which data file or output directories are used.

# Fields
- `data_dir`       : directory containing input TSV files
- `figures_dir`    : directory for generated plots (created automatically)
- `analysis_dir`   : directory for generated statistics tables (TSV)
- `metrics_file`   : path to the primary knot-seed-catalog file
- `cuts_file`      : path to the cuts-to-zero evaluation TSV
- `cut_eval_file`  : path to the single-cut evaluation TSV
- `cut_eval_stage_file`     : path to the single-cut stage metrics TSV
- `cuts_to_zero_stage_file` : path to the cuts-to-zero stage metrics TSV
- `complexity_dir`          : output directory for complexity outputs
- `complexity_features`     : metrics used to fit the complexity model
- `complexity_knot_filter`  : which knot types train the model
- `complexity_orientation`  : `:positive` (absolute PC1) or `:signed`
- `complexity_min_training_rows` : fall back to all rows below this sample size
- `correlation_metrics`     : metrics used in the correlation analysis

# Example
    s = Settings(metrics_file = "data/custom_catalog.tsv")
    run_all(s)
"""
Base.@kwdef struct Settings
    data_dir::String = joinpath(@__DIR__, "..", "data")
    figures_dir::String = joinpath(@__DIR__, "..", "figures")
    analysis_dir::String = joinpath(@__DIR__, "..", "analysis")
    metrics_file::String = joinpath(@__DIR__, "..", "data",
        "knot_seed_catalog_v3_1403_sorted_by_type_crossings.tsv")
    cuts_file::String = joinpath(@__DIR__, "..", "data", "knot_cuts_to_zero_v3.tsv")
    cut_eval_file::String = joinpath(@__DIR__, "..", "data", "knot_cut_eval_v3.tsv")
    cut_eval_stage_file::String = joinpath(@__DIR__, "..", "data",
        "knot_cut_eval_stage_metrics.tsv")
    cuts_to_zero_stage_file::String = joinpath(@__DIR__, "..", "data",
        "knot_cuts_to_zero_stage_metrics.tsv")
    complexity_dir::String = joinpath(@__DIR__, "..", "analysis", "complexity")
    complexity_features::Vector{String} = [
        "segmentCount", "selfCollisionCount", "closePairCount",
        "totalCurvature", "outsideKnotCurvature",
    ]
    complexity_knot_filter::String = "identified"
    complexity_orientation::Symbol = :positive
    complexity_min_training_rows::Int = 3
    correlation_metrics::Vector{String} = [
        "segmentCount", "totalHoopCount", "threadedHoopCount",
        "totalThreadingCount", "selfCollisionCount", "closePairCount",
        "outsideKnotCurvature",
    ]
end

"""
    DEFAULT_SETTINGS

The default configuration used by all functions when no `settings` argument
is provided.  Assign a new `Settings` instance to this variable to change
the defaults globally for the session:

    CableTangleAnalysis.DEFAULT_SETTINGS =
        CableTangleAnalysis.Settings(metrics_file = "other.tsv")
"""
const DEFAULT_SETTINGS = Settings()

# ---- Forward platform setup --------------------------------------------
using CairoMakie
CairoMakie.activate!(type = "pdf")
set_theme!(theme_minimal())

# ---- Load sub-module files (order matters) -----------------------------
include("data.jl")        # Data loading and knot grouping
include("utils.jl")       # TSV and column helpers
include("plotutils.jl")   # Shared plotting utilities
include("paired.jl")      # Paired statistics and violin plots
include("cuts.jl")        # Cut-evaluation analysis
include("complexity.jl")  # PCA-derived complexity score
include("pca.jl")         # PCA projection analysis
include("correlation.jl") # Correlation analysis
include("analysis.jl")    # Knot-type distribution

# ---- Exports -----------------------------------------------------------
export Settings, DEFAULT_SETTINGS, run_all,
       load_data, load_data_and_group_knots,
       clear_data_of_unwanted_knots,
       plot_cuts_histogram, plot_cut_eval_comparison,
       plot_threading_reduction_first_cut, plot_threading_reduction_resolved,
       plot_cuts_needed_all, plot_cuts_needed_resolved,
       first_cut_threading_reduction, resolved_threading_reduction,
       cuts_needed_summary,
       threading_reduction_first_cut_paired_tests,
       threading_reduction_resolved_paired_tests,
       cuts_needed_paired_tests, cuts_needed_resolved_paired_tests,
       cuts_violin_comparisons,
       ComplexityModel, ComplexityScores, COMPLEXITY_FEATURES,
       fit_complexity_pca, score_complexity, score_stage_metrics,
       analyze_complexity, plot_complexity_reduction_first_cut,
       plot_complexity_reduction_resolved,
       complexity_reduction_first_cut_paired_tests,
       complexity_reduction_resolved_paired_tests,
       resolved_complexity_reduction, stage_comparison_deltas,
       augment_table_with_complexity, write_complexity_model,
       write_complexity_scores, seed_score_map,
       run_pca, knot_distribution, correlation

# ---- Entry point -------------------------------------------------------

"""
    run_all(settings::Settings = DEFAULT_SETTINGS)

Run the complete analysis pipeline and write every figure and table to the
directories in `settings`.  Individual analyses can also be called
separately — see the function-level docstrings.
"""
function run_all(settings::Settings = DEFAULT_SETTINGS)
    println("═══ CableTangleAnalysis — Full Analysis Pipeline ═══\n")

    println("── 1. Cut evaluation ──")
    plot_cuts_histogram(settings.cuts_file; settings = settings)
    plot_cut_eval_comparison(settings.cut_eval_file; settings = settings)
    cuts_violin_comparisons(settings)

    println("\n── 2. Distribution, PCA and correlation ──")
    clear_data_of_unwanted_knots(settings)
    knot_distribution(settings)
    run_pca(settings)
    correlation(settings)

    println("\n── 3. Complexity analysis ──")
    model = analyze_complexity(settings)
    for (stage_file, plot_function, test_function) in (
        (settings.cut_eval_stage_file,
            plot_complexity_reduction_first_cut,
            complexity_reduction_first_cut_paired_tests),
        (settings.cuts_to_zero_stage_file,
            plot_complexity_reduction_resolved,
            complexity_reduction_resolved_paired_tests),
    )
        if isfile(stage_file)
            plot_function(stage_file; settings = settings, model = model)
            test_function(stage_file; settings = settings, model = model)
        else
            @warn "missing stage metrics; skipping complexity comparison" stage_file
        end
    end

    println("\n✅ All analyses complete. See $(settings.figures_dir)/ and $(settings.analysis_dir)/ for outputs.")
    return nothing
end

end # module CableTangleAnalysis
