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
    println("  julia mondrian.jl gen-prime [OPTIONS]       Generate A340984 prime rectangulations")
    println("  julia mondrian.jl -h, --help                Show this message")
    println()
    println("Options for gen-prime (run with `julia -t auto` for full throughput):")
    println("  --kmax <K>                 Highest piece count to grow to (default: 18)")
    println("  --per-k <NUM>              Distinct classes per generated k (default: 25000)")
    println("  --seconds-per-k <SEC>      Time budget per generated k (default: 900)")
    println("  --catalogue-k <K>          Highest complete catalogue level used (default: 14)")
    println("  --out <file>               Output library (default: data/import/prime_generated.cache)")
    println("  --jsonl <dir>              Also export A340984-schema .jsonl files to <dir>")
    println("  --prime-dir <dir>          A340984 representatives dir (default: data/import/a340984)")
    println()
    println("Options for init-cache:")
    println("  --expand-electrical <NUM>  Target electrical topologies per k (default: 20000)")
    println("  --expand-quilts <NUM>      Target quilt child topologies (default: 50000)")
    println("  --expand-prime <NUM>       Target A340984 prime topologies per k (default: 8000)")
    println("  --prime-per-k <NUM>        A340984 seed classes loaded per k, 0 = all (default: 8000)")
    println("  --prime-dir <dir>          A340984 representatives dir (default: <import>/a340984)")
    println("  --prime-generated <file>   Library from `gen-prime` (default: <import>/prime_generated.cache)")
    println("  --no-prime                 Skip the A340984 prime seeds entirely")
    println("  --fast-cache               Generate smaller cache (~20k topologies) for fast testing")
    println("  --cache-dir <dir>          Cache directory path (default: data/cache)")
    println("  --import-dir <dir>         Seed import directory (default: data/import)")
    println()
    println("Options for solve:")
    println("  --timeout <sec>            Total time limit in seconds (default: none / unlimited)")
    println("  --deep                     Perform deep 5-minute/config hierarchical search")
    println("  --output <file>            Output JSON path (default: solutions.json)")
end

# ─── Prime (A340984) Rectangulation Generation ────────────────────────────────

function counts_by_k(fps::Vector{Floorplan})::String
    by_k = Dict{Int, Int}()
    for fp in fps
        by_k[length(fp)] = get(by_k, length(fp), 0) + 1
    end
    return join(["k=$k:$(by_k[k])" for k in sort(collect(keys(by_k)))], ", ")
end

function run_gen_prime(; prime_dir::String      = joinpath("data", "import", "a340984"),
                         out_path::String       = joinpath("data", "import", "prime_generated.cache"),
                         kmax::Int              = 18,
                         per_k::Int             = 25000,
                         seconds_per_k::Float64 = 900.0,
                         catalogue_k::Int       = 14,
                         jsonl_dir::Union{String, Nothing} = nothing,
                         seed::Int              = 340984)
    isdir(prime_dir) || error("A340984 representatives not found at $prime_dir")
    rng = Xoshiro(seed)

    println("==================================================")
    println("  A340984 Prime Rectangulation Generator")
    println("==================================================")
    println("  Catalogue directory  : ", prime_dir)
    println("  Complete levels      : k <= ", catalogue_k)
    println("  Generated levels     : k = ", catalogue_k + 1, " .. ", kmax)
    println("  Target per level     : ", per_k, " classes  (", seconds_per_k, "s budget)")
    println("  Worker threads       : ", Threads.nthreads())
    Threads.nthreads() == 1 &&
        println("  [!] Single-threaded. Re-run as `julia -t auto mondrian.jl gen-prime` for ~$(Sys.CPU_THREADS)x throughput.")
    println()

    t0 = time()
    print("  [*] Loading the complete catalogue... ")
    base = parse_a340984(prime_dir; k_range=2:catalogue_k, max_per_k=0, rng=rng)
    @printf("%d classes in %.1fs\n", length(base), time() - t0)
    println("      ", counts_by_k(base))

    println("  [*] Growing new levels by split-and-repair...")
    on_step = (k, parents, found, elapsed) -> begin
        @printf("      [k=%-2d] %7d classes from %7d parents  (%.1fs, %.0f/s)\n",
                k, found, parents, elapsed, found / max(elapsed, 1e-9))
    end
    levels = generate_prime_ladder(base, kmax, per_k; rng=rng, seconds_per_k=seconds_per_k, on_step=on_step)

    generated = Floorplan[]
    for k in sort(collect(keys(levels)))
        append!(generated, levels[k])
    end

    if isempty(generated)
        println("  [!] No new levels generated.")
        return 0
    end

    sample = generated[1:max(1, cld(length(generated), 2000)):end]
    bad = count(fp -> !(is_valid(fp) && is_prime_rectangulation(fp)), sample)
    @printf("  [+] Audited %d of %d generated tilings: %d not prime/valid.\n", length(sample), length(generated), bad)

    save_library(out_path, generated)
    @printf("  [+] Saved %d generated classes to %s\n", length(generated), out_path)
    if jsonl_dir !== nothing
        files = save_a340984_jsonl(jsonl_dir, generated)
        @printf("  [+] Exported %d A340984-schema .jsonl files to %s\n", length(files), jsonl_dir)
    end
    println("      ", counts_by_k(generated))
    @printf("\nTotal prime library: %d classes (%d complete + %d generated) in %.1fs.\n",
            length(base) + length(generated), length(base), length(generated), time() - t0)
    return 0
end

# ─── Full-Scale Topology Cache Initialization ─────────────────────────────────

function run_init_cache(; import_dir::String        = joinpath("data", "import"),
                          cache_dir::String         = joinpath("data", "cache"),
                          expand_electrical::Int    = 20000,
                          expand_quilts::Int        = 50000,
                          expand_prime::Int         = 8000,
                          prime_per_k::Int          = 8000,
                          prime_dir::Union{String, Nothing} = nothing,
                          generated_path::String    = joinpath("data", "import", "prime_generated.cache"),
                          use_prime::Bool           = true,
                          prime_roots_per_r::Int    = 6,
                          seed::Int                 = 202609)
    rng = Xoshiro(seed)
    prime_path = prime_dir === nothing ? joinpath(import_dir, "a340984") : prime_dir
    prime_per_k <= 0 && (prime_per_k = typemax(Int))

    println("==================================================")
    println("  Full-Scale Mondrian Topology Cache Generator")
    println("==================================================")
    println("  Import directory     : ", import_dir)
    println("  Cache directory      : ", cache_dir)
    println("  Electrical expansion : ", expand_electrical, " per k")
    println("  Quilt expansion      : ", expand_quilts, " topologies")
    println("  Prime expansion      : ", use_prime ? "$(expand_prime) per k (<= $(prime_per_k) seeds/k)" : "disabled")
    println()

    known_path = joinpath(import_dir, "ed_pegg_jr.json")
    small_path = joinpath(import_dir, "oeis_solutions.txt")

    known_fps = isfile(known_path) ? import_floorplans(known_path) : Floorplan[]
    small_fps = isfile(small_path) ? parse_oeis_quilts(small_path; max_n=96) : Floorplan[]

    println("  [+] Ingested seeds: $(length(known_fps)) from ed_pegg_jr.json, $(length(small_fps)) from oeis_solutions.txt")

    # A340984: prime rectangulations, i.e. no proper sub-collection of 2..k-1
    # pieces forms a rectangle (so in particular no two pieces share a face).
    prime_fps = Floorplan[]
    if use_prime
        if isdir(prime_path)
            prime_fps = parse_a340984(prime_path; k_range=5:18, max_per_k=prime_per_k, rng=rng)
        else
            println("  [!] No A340984 representatives at $(prime_path) - catalogue seeds skipped.")
        end

        # Levels beyond the catalogue, if `mondrian.jl gen-prime` has been run.
        if isfile(generated_path)
            pooled = Dict{Int, Vector{Floorplan}}()
            for fp in load_library(generated_path)
                bucket = get!(pooled, length(fp), Floorplan[])
                length(bucket) < prime_per_k && push!(bucket, fp)
            end
            extra = vcat(values(pooled)...)
            println("  [+] Loaded $(length(extra)) generated prime classes from $(generated_path)")
            append!(prime_fps, extra)
            sort!(prime_fps, by=fp -> (length(fp), fp.id))
        end
    end
    isempty(prime_fps) ||
        println("  [+] A340984 prime seeds: $(length(prime_fps)) classes  ($(counts_by_k(prime_fps)))")

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

    expanded_prime = Floorplan[]
    if !isempty(prime_fps)
        println("  [*] Expanding A340984 prime seed caches (target: $(expand_prime) per k)...")
        on_step_prime = (k, initial_n, final_n, elapsed) -> begin
            @printf("      [k=%-2d] %4d -> %-5d topologies  (%.2fs)\n", k, initial_n, final_n, elapsed)
        end
        expanded_prime = generate_prime_portfolio(prime_fps, expand_prime; k_range=8:18, rng=rng, on_step=on_step_prime)
        println("  [+] Total expanded prime topologies: $(length(expanded_prime)).")
    end

    # The small prime tilings (pinwheel and friends) enter the quilt engine as
    # substitution blocks: they are irreducible, so they can only ever appear in
    # a larger rectangulation by being substituted in.
    prime_blocks = Floorplan[]
    for fp in prime_fps
        2 <= length(fp) <= 8 || continue
        eq = equalize(fp)
        eq !== nothing && push!(prime_blocks, eq)
    end

    println("  [*] Generating multi-level quilt portfolio from oeis_solutions.txt (target: $(expand_quilts))...")
    isempty(prime_blocks) || println("      (+ $(length(prime_blocks)) prime substitution blocks)")
    on_step_quilt = (depth, curr_count, target, elapsed) -> begin
        @printf("      [depth=%-2d] %5d / %-5d topologies  (%.2fs)\n", depth, curr_count, target, elapsed)
    end
    quilt_portfolio = generate_quilt_portfolio(small_fps, expand_quilts; mink=10, maxk=18,
                                               extra_blocks=prime_blocks, rng=rng, on_step=on_step_quilt)
    println("  [+] Total quilt child portfolio: $(length(quilt_portfolio)) topologies.")

    println("  [*] Building outer root templates (R = 2..14)...")
    outer_roots = Floorplan[]
    push!(outer_roots, Floorplan("outer_root_R2_1", [(0.0, 0.0, 0.5, 1.0), (0.5, 0.0, 1.0, 1.0)]))

    known_by_r = Dict{Int, Vector{Floorplan}}()
    for fp in vcat(known_fps, small_fps)
        4 <= length(fp) <= 14 || continue
        eq = equalize(fp)
        eq === nothing && continue
        push!(get!(known_by_r, length(fp), Floorplan[]), Floorplan("outer_known_$(fp.id)", eq.rectangles))
    end

    prime_by_r = Dict{Int, Vector{Floorplan}}()
    for fp in prime_fps
        4 <= length(fp) <= 14 || continue
        length(get(prime_by_r, length(fp), Floorplan[])) >= prime_roots_per_r && continue
        eq = equalize(fp)
        eq === nothing && continue
        push!(get!(prime_by_r, length(fp), Floorplan[]), Floorplan("outer_$(fp.id)", eq.rectangles))
    end

    # solve() only tries the first few roots of each R, so the two families are
    # interleaved per R instead of concatenated.
    known_roots_count = sum(length, values(known_by_r); init=0)
    prime_roots_count = sum(length, values(prime_by_r); init=0)
    for r in 2:14
        append!(outer_roots, interleave(get(known_by_r, r, Floorplan[]),
                                        get(prime_by_r, r, Floorplan[])))
    end

    for r in 3:8
        raw_roots = generate_slicing_library(12; k_range=r:r, rng=rng)
        for (idx, fp) in enumerate(raw_roots)
            eq = equalize(fp)
            eq !== nothing && push!(outer_roots, Floorplan("outer_slicing_R$(r)_$(idx)", eq.rectangles))
        end
    end
    println("  [+] Total equalized outer root templates: $(length(outer_roots)) ($(known_roots_count) known, $(prime_roots_count) prime).")

    mkpath(cache_dir)

    # Interleave the engines within each piece count: the single-level search
    # walks children in stored order under a time cap, so a plain vcat would
    # bury whichever engine comes last.
    full_child_library = Floorplan[]
    sources = [expanded_electrical, quilt_portfolio, expanded_prime]
    for k in sort(collect(Set(length(fp) for src in sources for fp in src)))
        append!(full_child_library, interleave([[fp for fp in src if length(fp) == k] for src in sources]...))
    end
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

    if args[1] == "gen-prime"
        prime_dir     = joinpath("data", "import", "a340984")
        out_path      = joinpath("data", "import", "prime_generated.cache")
        kmax          = 18
        per_k         = 25000
        seconds_per_k = 900.0
        catalogue_k   = 14
        jsonl_dir     = nothing

        i = 2
        while i <= length(args)
            if args[i] == "--prime-dir"; i += 1; prime_dir = args[i]
            elseif args[i] == "--out"; i += 1; out_path = args[i]
            elseif args[i] == "--kmax"; i += 1; kmax = parse(Int, args[i])
            elseif args[i] == "--per-k"; i += 1; per_k = parse(Int, args[i])
            elseif args[i] == "--seconds-per-k"; i += 1; seconds_per_k = parse(Float64, args[i])
            elseif args[i] == "--catalogue-k"; i += 1; catalogue_k = parse(Int, args[i])
            elseif args[i] == "--jsonl"; i += 1; jsonl_dir = args[i]
            else
                println(stderr, "Unknown argument: ", args[i]); return 1
            end
            i += 1
        end
        return run_gen_prime(; prime_dir=prime_dir, out_path=out_path, kmax=kmax, per_k=per_k,
                               seconds_per_k=seconds_per_k, catalogue_k=catalogue_k, jsonl_dir=jsonl_dir)
    end

    if args[1] == "init-cache"
        import_dir         = joinpath("data", "import")
        cache_dir          = joinpath("data", "cache")
        expand_electrical  = 20000
        expand_quilts      = 50000
        expand_prime       = 8000
        prime_per_k        = 8000
        prime_dir          = nothing
        generated_path     = joinpath("data", "import", "prime_generated.cache")
        use_prime          = true

        i = 2
        while i <= length(args)
            if args[i] == "--import-dir"; i += 1; import_dir = args[i]
            elseif args[i] == "--cache-dir"; i += 1; cache_dir = args[i]
            elseif args[i] == "--expand-electrical"; i += 1; expand_electrical = parse(Int, args[i])
            elseif args[i] == "--expand-quilts"; i += 1; expand_quilts = parse(Int, args[i])
            elseif args[i] == "--expand-prime"; i += 1; expand_prime = parse(Int, args[i])
            elseif args[i] == "--prime-per-k"; i += 1; prime_per_k = parse(Int, args[i])
            elseif args[i] == "--prime-dir"; i += 1; prime_dir = args[i]
            elseif args[i] == "--prime-generated"; i += 1; generated_path = args[i]
            elseif args[i] == "--no-prime"; use_prime = false
            elseif args[i] == "--fast-cache"
                expand_electrical = 2500
                expand_quilts     = 8000
                expand_prime      = 1500
                prime_per_k       = 300
            else
                println(stderr, "Unknown argument: ", args[i]); return 1
            end
            i += 1
        end
        run_init_cache(; import_dir=import_dir, cache_dir=cache_dir,
                         expand_electrical=expand_electrical, expand_quilts=expand_quilts,
                         expand_prime=expand_prime, prime_per_k=prime_per_k,
                         prime_dir=prime_dir, generated_path=generated_path, use_prime=use_prime)
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