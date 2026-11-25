#!/usr/bin/env julia
#
# Run the CableTangleAnalysis pipeline for a CableSim dataset.
#
# Usage:
#   julia --project=<repo>/cabletangleanalysis \
#       <repo>/cabletangleanalysis/scripts/run_dataset.jl <dataset-root> [variant ...]
#
# A variant is either a name below the dataset root (`zero`), a directory
# (`path/to/variant`), or `label=path` when the label should differ from the
# directory name.  For a flat dataset whose tables sit at the root:
#   run_dataset.jl results/dataset_final single=results/dataset_final
#
# Each variant directory must contain catalog.tsv, cut_eval.tsv,
# cuts_to_zero.tsv, cut_eval_stage_metrics.tsv, and
# cuts_to_zero_stage_metrics.tsv.  Outputs are written to
# <dataset-root>/analysis_julia/<label>/.

using CableTangleAnalysis
const CTA = CableTangleAnalysis

const ORIENTATION = Symbol(get(ENV, "CABLESIM_PCA_ORIENTATION", "positive"))
const COMPLEXITY_FILTER = get(ENV, "CABLESIM_COMPLEXITY_FILTER", "identified")
const PER_VARIANT_COMPLEXITY = get(ENV, "CABLESIM_COMPLEXITY_PER_VARIANT", "1") != "0"

function resolve_variant(dataset_root::String, variant::String)
    # Accept a label, a directory path, or "label=/path/to/variant".
    if occursin("=", variant)
        parts = split(variant, "=", limit = 2)
        if isdir(parts[2])
            label = isempty(parts[1]) ? basename(parts[2]) : String(parts[1])
            return label, String(parts[2])
        end
    end
    isdir(variant) && return basename(variant), variant
    return variant, joinpath(dataset_root, variant)
end

function variant_settings(dataset_root::String, variant::String)
    label, variant_dir = resolve_variant(dataset_root, variant)
    isdir(variant_dir) || error(
        "variant directory not found: $variant_dir\n" *
        "usage: run_dataset.jl <dataset-root> [variant ...]\n" *
        "for a flat dataset use: run_dataset.jl <dataset-root> single=<dataset-root>",
    )
    output_dir = joinpath(dataset_root, "analysis_julia", label)
    settings = CTA.Settings(
        data_dir = joinpath(output_dir, "data"),
        figures_dir = joinpath(output_dir, "figures"),
        analysis_dir = joinpath(output_dir, "tables"),
        metrics_file = joinpath(variant_dir, "catalog.tsv"),
        cuts_file = joinpath(variant_dir, "cuts_to_zero.tsv"),
        cut_eval_file = joinpath(variant_dir, "cut_eval.tsv"),
        cut_eval_stage_file = joinpath(variant_dir, "cut_eval_stage_metrics.tsv"),
        cuts_to_zero_stage_file = joinpath(variant_dir, "cuts_to_zero_stage_metrics.tsv"),
        complexity_dir = joinpath(output_dir, "complexity"),
        complexity_orientation = ORIENTATION,
        complexity_knot_filter = COMPLEXITY_FILTER,
    )
    mkpath(settings.data_dir)
    return label, settings
end

function run_variant(dataset_root::String, variant::String)
    label, settings = variant_settings(dataset_root, variant)

    println("\n=== Cut evaluation for $label ===")
    CTA.plot_cuts_histogram(settings.cuts_file; settings = settings)
    CTA.plot_cut_eval_comparison(settings.cut_eval_file; settings = settings)
    CTA.cuts_violin_comparisons(settings)

    println("\n=== Distribution, PCA and correlation for $label ===")
    CTA.clear_data_of_unwanted_knots(settings)
    CTA.knot_distribution(settings)
    CTA.run_pca(settings)
    CTA.correlation(settings)

    PER_VARIANT_COMPLEXITY || return nothing

    println("\n=== Complexity analysis for $label ===")
    model = CTA.analyze_complexity(settings)
    for (stage_file, plot_function, test_function) in (
        (settings.cut_eval_stage_file,
            CTA.plot_complexity_reduction_first_cut,
            CTA.complexity_reduction_first_cut_paired_tests),
        (settings.cuts_to_zero_stage_file,
            CTA.plot_complexity_reduction_resolved,
            CTA.complexity_reduction_resolved_paired_tests),
    )
        if isfile(stage_file)
            plot_function(stage_file; settings = settings, model = model)
            test_function(stage_file; settings = settings, model = model)
        else
            @warn "missing stage metrics; skipping complexity comparison" stage_file
        end
    end
    return nothing
end

# Pooled complexity model across variants.  The per-variant models above stay
# independent; this one exists for cross-variant comparisons.
function run_pooled(dataset_root::String, variants::Vector{String})
    catalogs = String[]
    for variant in variants
        _, variant_dir = resolve_variant(dataset_root, variant)
        catalog = joinpath(variant_dir, "catalog.tsv")
        isfile(catalog) && push!(catalogs, catalog)
    end
    length(catalogs) >= 2 || return nothing

    println("\n=== Pooled complexity model ===")
    output_dir = joinpath(dataset_root, "complexity")
    model, records = CTA.fit_complexity_pca(
        catalogs;
        features = CTA.COMPLEXITY_FEATURES,
        knot_filter = COMPLEXITY_FILTER,
        orientation = ORIENTATION,
    )
    CTA.write_complexity_model(model, output_dir)
    CTA.write_complexity_scores(model, records, output_dir)
    println("Pooled model written to $output_dir")
    return nothing
end

function main(args)
    isempty(args) && error("usage: run_dataset.jl <dataset-root> [variant ...]")
    first_arg = args[1]
    if occursin("=", first_arg)
        error(
            "usage: run_dataset.jl <dataset-root> [variant ...]\n" *
            "`$first_arg` is a variant spec, not a dataset root. " *
            "Pass the root first, for example:\n" *
            "  run_dataset.jl <dataset-root> $first_arg",
        )
    end
    dataset_root = abspath(first_arg)
    isdir(dataset_root) || error("dataset root is not a directory: $dataset_root")
    variants = length(args) > 1 ? String.(args[2:end]) : ["aggregate"]

    for variant in variants
        run_variant(dataset_root, variant)
    end
    run_pooled(dataset_root, variants)
    println("\nAll requested variants complete.")
end

main(ARGS)
