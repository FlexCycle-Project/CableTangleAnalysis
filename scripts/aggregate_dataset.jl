#!/usr/bin/env julia
#
# Aggregate per-seed result directories into the aggregate TSVs expected by
# cabletangleanalysis/scripts/run_dataset.jl, then optionally run the analysis.
#
# Usage:
#   julia --project=cabletangleanalysis scripts/aggregate_dataset.jl <dataset-root> [label] [--no-analysis]
#
# Layout:
#   <root>/seed_<N>/{catalog.tsv,cut_eval.tsv,cuts_to_zero.tsv,
#                    cut_eval_stage_metrics.tsv,cuts_to_zero_stage_metrics.tsv,
#                    topology_audit_*.tsv}
#   <root>/<label>/...          aggregate tables written here
#   <root>/analysis_julia/<label>/...  analysis output

struct TableSpec
    name::String
    required::Bool
end

function usage()
    println(stderr, "usage: aggregate_dataset.jl <dataset-root> [label] [--no-analysis]")
    exit(2)
end

function first_line(path::String)
    open(path, "r") do io
        return readline(io)
    end
end

function append_rows!(out_io, path::String, expected_header::String)
    open(path, "r") do in
        header = readline(in)
        header == expected_header || error("header mismatch in $path")
        for line in eachline(in)
            isempty(line) || println(out_io, line)
        end
    end
end

function read_meta(path::String)
    meta = Dict{String, String}()
    for line in eachline(path)
        isempty(strip(line)) && continue
        parts = split(line, '\t'; limit = 2)
        length(parts) == 2 && (meta[parts[1]] = parts[2])
    end
    return meta
end

# Every seed must come from the same generator configuration; otherwise the
# aggregate mixes incomparable runs.
function validate_metas(seed_dirs)
    config_keys = [
        "version", "creep", "creep_settle_steps", "creep_threshold", "dt",
        "random_trials", "ctz_max_cuts", "ctz_random_trials",
        "audit_steps", "audit_cuts", "audit_contact_eps", "audit_trace",
        "oracle", "audit",
    ]
    metas = Dict{String, Dict{String, String}}()
    complete_dirs = String[]
    for dir in seed_dirs
        path = joinpath(dir, "meta.tsv")
        if !isfile(path)
            @warn "skipping seed directory without meta.tsv: $dir"
            continue
        end
        meta = read_meta(path)
        if get(meta, "status", "") != "complete"
            @warn "skipping incomplete seed $(basename(dir)) (status=$(get(meta, "status", "missing")))"
            continue
        end
        push!(complete_dirs, dir)
        metas[basename(dir)] = meta
    end
    isempty(metas) && error("no complete seed runs found")

    reference_name = first(sort(collect(keys(metas))))
    reference = metas[reference_name]
    for key in config_keys
        expected = get(reference, key, "")
        for (name, meta) in metas
            actual = get(meta, key, "")
            actual == expected || error(
                "configuration mismatch for '$key': $reference_name=$expected, $name=$actual",
            )
        end
    end
    return complete_dirs
end

function main(args)
    root = ""
    label = "aggregate"
    run_analysis = true

    for arg in args
        if arg == "--no-analysis"
            run_analysis = false
        elseif startswith(arg, "-")
            usage()
        elseif isempty(root)
            root = arg
        elseif label == "aggregate"
            label = arg
        else
            usage()
        end
    end
    isempty(root) && usage()
    root = abspath(root)
    isdir(root) || error("dataset root is not a directory: $root")

    outdir = joinpath(root, label)
    mkpath(outdir)

    seed_dirs = sort(
        filter(
            d -> isdir(d) && startswith(basename(d), "seed_"),
            readdir(root, join = true),
        ),
        by = d -> try
            parse(Int, split(basename(d), "_"; limit = 2)[2])
        catch
            typemax(Int)
        end,
    )
    isempty(seed_dirs) && error("no seed_* directories under $root")
    seed_dirs = validate_metas(seed_dirs)

    tables = [
        TableSpec("catalog.tsv", true),
        TableSpec("cut_eval.tsv", true),
        TableSpec("cuts_to_zero.tsv", true),
        TableSpec("cut_eval_stage_metrics.tsv", true),
        TableSpec("cuts_to_zero_stage_metrics.tsv", true),
        TableSpec("topology_audit_seeds.tsv", false),
        TableSpec("topology_audit_events.tsv", false),
    ]

    for table in tables
        paths = [
            joinpath(d, table.name) for d in seed_dirs
            if isfile(joinpath(d, table.name))
        ]
        if isempty(paths)
            table.required && error(
                "missing required table $(table.name) in all seed directories",
            )
            println("[aggregate] skip $(table.name) (absent)")
            continue
        end
        if length(paths) != length(seed_dirs)
            @warn "only $(length(paths))/$(length(seed_dirs)) seeds have $(table.name)"
        end

        header = first_line(paths[1])
        for path in paths[2:end]
            first_line(path) == header || error(
                "header mismatch between $(paths[1]) and $path",
            )
        end

        destination = joinpath(outdir, table.name)
        open(destination, "w") do out
            println(out, header)
            for path in paths
                append_rows!(out, path, header)
            end
        end
        println("[aggregate] $(table.name): $(length(paths)) seed file(s) -> $destination")
    end

    println("[aggregate] aggregate tables written to $outdir")

    if run_analysis
        script_dir = @__DIR__
        project = normpath(joinpath(script_dir, "..", "cabletangleanalysis"))
        run_dataset = joinpath(project, "scripts", "run_dataset.jl")
        isfile(run_dataset) || error("analysis driver not found: $run_dataset")

        cmd = `$(Base.julia_cmd()) --project=$project $run_dataset $root $label=$outdir`
        println("[aggregate] running: ", cmd)
        run(cmd)
    end
end

main(ARGS)
