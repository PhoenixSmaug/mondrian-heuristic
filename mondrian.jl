#!/usr/bin/env julia

# ==============================================================================
# mondrian.jl — Master Mondrian Puzzle Solver & Cache CLI
# ==============================================================================

if !isdefined(@__MODULE__, :MondrianSolver)
    include(joinpath(@__DIR__, "src", "MondrianSolver.jl"))
end

using .MondrianSolver
using Random
using Printf

function print_usage()
    println("Mondrian Art Problem Solver CLI")
    println()
    println("Usage:")
    println("  julia mondrian.jl <N> [OPTIONS]             Solve for a single grid size N")
    println("  julia mondrian.jl <N_START> <N_END> [OPT]   Solve for an inclusive range of N")
    println("  julia mondrian.jl init-cache [OPTIONS]      Generate full topology cache")
    println("  julia mondrian.jl -h, --help                Show this message")
    println()
    println("Options for init-cache:")
    println("  --expand-electrical <NUM>  Target electrical topologies per k (default: 20000)")
    println("  --expand-quilts <NUM>      Target quilt child topologies (default: 50000)")
    println("  --fast-cache               Generate smaller cache (~20k topologies) for fast testing")
    println("  --cache-dir <dir>          Cache directory path (default: data/cache)")
    println("  --import-dir <dir>         Seed import directory (default: data/import)")
    println()
    println("Options for solve:")
    println("  --timeout <sec>            Total time limit in seconds (default: none / unlimited)")
    println("  --deep                     Perform deep 5-minute/config hierarchical search")
    println("  --output <file>            Output JSON path (default: solutions.json)")
end

# ─── Full-Scale Topology Cache Initialization ─────────────────────────────────

function run_init_cache(; import_dir::String        = joinpath("data", "import"),
                          cache_dir::String         = joinpath("data", "cache"),
                          expand_electrical::Int    = 20000,
                          expand_quilts::Int        = 50000,
                          seed::Int                 = 202609)
    rng = Xoshiro(seed)
    println("==================================================")
    println("  Full-Scale Mondrian Topology Cache Generator")
    println("==================================================")
    println("  Import directory     : ", import_dir)
    println("  Cache directory      : ", cache_dir)
    println("  Electrical expansion : ", expand_electrical, " per k")
    println("  Quilt expansion      : ", expand_quilts, " topologies")
    println()

    known_path = joinpath(import_dir, "ed_pegg_jr.json")
    small_path = joinpath(import_dir, "oeis_solutions.txt")

    known_fps = isfile(known_path) ? import_floorplans(known_path) : Floorplan[]
    small_fps = isfile(small_path) ? parse_oeis_quilts(small_path; max_n=96) : Floorplan[]

    println("  [+] Ingested seeds: $(length(known_fps)) from ed_pegg_jr.json, $(length(small_fps)) from oeis_solutions.txt")

    println("  [*] Synthesizing base planar electrical dual seeds (k=10..15)...")
    base_electrical = Dict{Int, Vector{Floorplan}}()
    for k in 10:15
        base_electrical[k] = generate_electrical_library(150; k_range=k:k, rng=rng)
        @printf("      [k=%-2d] %4d base electrical seeds synthesized\n", k, length(base_electrical[k]))
    end

    println("  [*] Expanding electrical dual caches (target: $(expand_electrical) per k)...")
    on_step_elec = (k, initial_n, final_n, elapsed) -> begin
        @printf("      [k=%-2d] %4d -> %-5d topologies  (%.2fs)\n", k, initial_n, final_n, elapsed)
    end
    expanded_electrical = expand_electrical_portfolio(base_electrical, expand_electrical; rng=rng, on_step=on_step_elec)
    println("  [+] Total expanded electrical duals: $(length(expanded_electrical)) topologies.")

    println("  [*] Generating multi-level quilt portfolio from oeis_solutions.txt (target: $(expand_quilts))...")
    on_step_quilt = (depth, curr_count, target, elapsed) -> begin
        @printf("      [depth=%-2d] %5d / %-5d topologies  (%.2fs)\n", depth, curr_count, target, elapsed)
    end
    quilt_portfolio = generate_quilt_portfolio(small_fps, expand_quilts; mink=10, maxk=18, rng=rng, on_step=on_step_quilt)
    println("  [+] Total quilt child portfolio: $(length(quilt_portfolio)) topologies.")

    println("  [*] Building outer root templates (R = 2..14)...")
    outer_roots = Floorplan[]
    push!(outer_roots, Floorplan("outer_root_R2_1", [(0.0, 0.0, 0.5, 1.0), (0.5, 0.0, 1.0, 1.0)]))

    known_roots_count = 0
    all_known_and_small = vcat(known_fps, small_fps)
    for fp in all_known_and_small
        k_fp = length(fp)
        if 4 <= k_fp <= 14
            eq = equalize(fp)
            if eq !== nothing
                push!(outer_roots, Floorplan("outer_known_$(fp.id)", eq.rectangles))
                known_roots_count += 1
            end
        end
    end

    for r in 3:8
        raw_roots = generate_slicing_library(12; k_range=r:r, rng=rng)
        for (idx, fp) in enumerate(raw_roots)
            eq = equalize(fp)
            eq !== nothing && push!(outer_roots, Floorplan("outer_slicing_R$(r)_$(idx)", eq.rectangles))
        end
    end
    println("  [+] Total equalized outer root templates: $(length(outer_roots)) ($(known_roots_count) non-slicing seeds).")

    mkpath(cache_dir)
    full_child_library = vcat(expanded_electrical, quilt_portfolio)
    for k in 8:16
        append!(full_child_library, generate_slicing_library(50; k_range=k:k, rng=rng))
    end

    lib_cache  = joinpath(cache_dir, "floorplans.cache")
    lib_json   = joinpath(cache_dir, "floorplans.json")
    root_cache = joinpath(cache_dir, "outer_roots.cache")

    print("  [*] Saving $(length(full_child_library)) child topologies to $(lib_cache)... ")
    save_library(lib_cache, full_child_library)
    println("done.")

    print("  [*] Saving $(length(outer_roots)) outer roots to $(root_cache)... ")
    save_library(root_cache, outer_roots)
    println("done.")

    println("\nCache generation complete: $(length(full_child_library)) child topologies, $(length(outer_roots)) outer templates.")
end

# ─── Solver Runner ────────────────────────────────────────────────────────────

function run_solve(targets::Vector{Int};
                   timeout::Float64                    = Inf,
                   deep::Bool                          = false,
                   output_path::Union{String, Nothing} = nothing,
                   cache_dir::String                   = joinpath("data", "cache"))

    lib_path  = joinpath(cache_dir, "floorplans.cache")
    root_path = joinpath(cache_dir, "outer_roots.cache")

    library = isfile(lib_path) ? load_library(lib_path) : Floorplan[]
    roots   = isfile(root_path) ? load_library(root_path) : Floorplan[]

    if isempty(library)
        println("Notice: Caches not found at $cache_dir. Run 'julia mondrian.jl init-cache' first.\n")
    end

    out_file = output_path === nothing ? "solutions.json" : output_path
    timeout_str = isfinite(timeout) ? "$(timeout)s" : "none"
    depth_str   = deep ? "DEEP (5 min/config)" : "standard (1 min/config)"

    println("==================================================")
    println("  Mondrian Art Problem Solver")
    println("  Target(s)    : $(targets[1])$(length(targets) > 1 ? " .. $(targets[end])" : "")")
    println("  Mode         : $depth_str")
    println("  Timeout      : $timeout_str")
    println("  Output File  : $out_file (consolidated, keeps best)")
    println("==================================================\n")

    start_t = time()
    is_interactive = isa(stdout, Base.TTY)

    for n in targets
        @printf("─── Solving N = %d ───\n", n)

        on_config_start = (c_idx, total_c, r, k, pcs, best_str, max_sec) -> begin
            is_interactive && print("\r" * " "^85 * "\r")
            @printf("  [*] [%5.1fs] Config %d/%d: R=%d, M=%d (Target pieces: %d | Incumbent: %s | Budget: %ds)\n",
                    time() - start_t, c_idx, total_c, r, k, pcs, best_str, max_sec)
            flush(stdout)
        end

        on_step_progress = (step, total_steps, root_i, total_roots, trial_i, total_trials, elapsed_cfg, max_sec, cap) -> begin
            if is_interactive
                pct = clamp(round(Int, (step / max(1, total_steps)) * 100), 0, 100)
                filled = div(pct, 5)
                bar = string("[", repeat("=", filled), (filled < 20 ? ">" : ""), repeat(" ", max(0, 20 - filled - (filled < 20 ? 1 : 0))), "]")
                cap_str = cap == typemax(Int) ? "none" : string(cap)
                @printf("\r      %s %3d%%  Root %d/%d | Trial %d/%d | %4.1fs/%-3ds | Cap: %-4s",
                        bar, pct, root_i, total_roots, trial_i, total_trials, elapsed_cfg, max_sec, cap_str)
                flush(stdout)
            end
        end

        on_update = (sol) -> begin
            is_interactive && print("\r" * " "^85 * "\r")
            updated, _ = update_solutions_json(out_file, sol)
            if updated
                @printf("  [+] [%5.1fs] NEW BEST defect = %-3d  (pieces: %-2d, %s) -> updated %s\n",
                        time() - start_t, sol.defect, length(sol.rectangles), sol.method, out_file)
            end
            flush(stdout)
        end

        sol = try
            solve(n; library=library, outer_roots=roots, timeout=timeout, deep=deep,
                  on_update=on_update, on_config_start=on_config_start, on_step_progress=on_step_progress)
        catch e
            if e isa InterruptException
                is_interactive && println()
                println("\n[!] Search interrupted by user.")
                println("[!] Current best solutions preserved in: ", out_file)
                return 0
            else
                rethrow(e)
            end
        end

        is_interactive && print("\r" * " "^85 * "\r")

        if sol !== nothing
            update_solutions_json(out_file, sol)
            @printf("SUCCESS N=%d | defect=%d | pieces=%d | method=%s\n\n",
                    sol.n, sol.defect, length(sol.rectangles), sol.method)
        else
            @printf("FAILED N=%d within timeout limit.\n\n", n)
        end
    end

    println("==================================================")
    println("Done: All solutions saved to: $out_file")
    println("Run 'julia verify.jl $out_file' to audit results.")
    return 0
end

# ─── Main CLI Dispatcher ──────────────────────────────────────────────────────

function main(args::Vector{String}=ARGS)
    if isempty(args) || args[1] in ("-h", "--help")
        print_usage(); return 0
    end

    if args[1] == "init-cache"
        import_dir         = joinpath("data", "import")
        cache_dir          = joinpath("data", "cache")
        expand_electrical  = 20000
        expand_quilts      = 50000

        i = 2
        while i <= length(args)
            if args[i] == "--import-dir"; i += 1; import_dir = args[i]
            elseif args[i] == "--cache-dir"; i += 1; cache_dir = args[i]
            elseif args[i] == "--expand-electrical"; i += 1; expand_electrical = parse(Int, args[i])
            elseif args[i] == "--expand-quilts"; i += 1; expand_quilts = parse(Int, args[i])
            elseif args[i] == "--fast-cache"
                expand_electrical = 2500
                expand_quilts     = 8000
            else
                println(stderr, "Unknown argument: ", args[i]); return 1
            end
            i += 1
        end
        run_init_cache(; import_dir=import_dir, cache_dir=cache_dir,
                         expand_electrical=expand_electrical, expand_quilts=expand_quilts)
        return 0
    end

    targets     = Int[]
    timeout     = Inf
    deep        = false
    output_path = nothing
    cache_dir   = joinpath("data", "cache")

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("-h", "--help"); print_usage(); return 0
        elseif arg == "--timeout"; i += 1; timeout = parse(Float64, args[i])
        elseif arg == "--deep"; deep = true
        elseif arg == "--output"; i += 1; output_path = args[i]
        elseif arg == "--cache-dir"; i += 1; cache_dir = args[i]
        elseif tryparse(Int, arg) !== nothing; push!(targets, parse(Int, arg))
        else; println(stderr, "Unknown argument: ", arg); return 1
        end
        i += 1
    end

    isempty(targets) && (println(stderr, "Error: Specify grid size N."); return 1)

    target_list = length(targets) == 1 ? targets : collect(targets[1]:targets[2])
    return run_solve(target_list; timeout=timeout, deep=deep, output_path=output_path, cache_dir=cache_dir)
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    exit(main(ARGS))
end