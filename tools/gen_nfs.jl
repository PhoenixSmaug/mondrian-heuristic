#!/usr/bin/env julia
# EXPERIMENTAL: not used for any result in the README. Writes
# data/import/nfs_generated.cache, which nothing else reads.
#
# Generate NON-FACE-SHARING rectangulations and measure the class.
#
# A rectangulation is face-share-free when no two of its rectangles are flush
# along a complete shared edge, i.e. no PAIR has a rectangular union. That is
# strictly weaker than being an A340984 prime rectangulation, which forbids
# EVERY sub-collection of 2..k-1 from tiling a rectangle. Prime implies
# face-share-free; the converse fails often.
#
# Why the weaker class is the relevant one for Mondrian: a shared complete
# face locks two rectangles' side lengths to each other, so the pair can only
# trade area through the single wall between them. A larger block merely
# constrains its own outline and leaves its interior free. It is the pairwise
# rigidity that blocks area equalization.
#
# Measured on the 2,991 reference dissections: 302 face-share-free, of which
# only 126 prime -- so 58% of the class is unreachable through primality.
#
#   julia -t auto tools/gen_nfs.jl [--per-k N] [--kmax K] [--seconds S]
#
# This file deliberately does NOT modify the solver module: it includes it and
# defines its own generator, so a mistake here cannot break the evaluator.

include(joinpath(@__DIR__, "..", "src", "MondrianSolver.jl"))
using .MondrianSolver
using Random

"""True when two vertices share a complete face, sides of the square included.

The four sides are vertices of the adjacency graph, so a rectangle covering
an entire side shares a complete face with that side vertex exactly as two
flush rectangles share one. A rectangle does that iff it spans a full
dimension while touching a parallel side -- which is also the only way a
rectangle can touch three sides of the square.
"""
function has_face_share(fp::Floorplan)::Bool
    r = fp.rectangles
    @inbounds for a in r
        (a[2] == 0.0 && a[4] == 1.0 && (a[1] == 0.0 || a[3] == 1.0)) && return true
        (a[1] == 0.0 && a[3] == 1.0 && (a[2] == 0.0 || a[4] == 1.0)) && return true
    end
    @inbounds for i in 1:length(r), j in 1:i-1
        a, b = r[i], r[j]
        (a[1] == b[1] && a[3] == b[3] && (a[4] == b[2] || b[4] == a[2])) && return true
        (a[2] == b[2] && a[4] == b[4] && (a[3] == b[1] || b[3] == a[1])) && return true
    end
    return false
end

is_face_share_free(fp::Floorplan)::Bool = !has_face_share(fp)

"""Split one rectangle, then repair by T-flips until face-share-free."""
function grow_nfs(fp::Floorplan, rng::AbstractRNG; max_repair::Int=4)
    cand = split_rectangle(fp, rand(rng, 1:length(fp)), rand(rng, Bool),
                           0.25 + 0.5 * rand(rng))
    cur = equalize(cand)
    cur === nothing && return nothing
    for _ in 1:max_repair
        is_face_share_free(cur) && return Floorplan("nfs_k$(length(cur))", cur.rectangles)
        flips = t_flips(cur)
        isempty(flips) && return nothing
        nxt = equalize(flips[rand(rng, 1:length(flips))])
        nxt === nothing && return nothing
        cur = nxt
    end
    # Also test the state after the final flip, which the loop never checks.
    return is_face_share_free(cur) ? Floorplan("nfs_k$(length(cur))", cur.rectangles) : nothing
end

function generate_level(pool::Vector{Floorplan}, target::Int;
                        max_seconds::Float64=120.0, max_repair::Int=4)
    (isempty(pool) || target <= 0) && return Floorplan[], 0
    nthreads = max(1, Threads.nthreads())
    deadline = time() + max_seconds
    seeds = [rand(UInt64) for _ in 1:nthreads]
    per_worker = cld(11 * target, 10 * nthreads)
    parts = [Floorplan[] for _ in 1:nthreads]
    tries = zeros(Int, nthreads)
    Threads.@threads for w in 1:nthreads
        # Distinct NAMES matter here. Julia scopes a variable assigned inside
        # a loop to the enclosing function if that name is also assigned at
        # function level -- and `seen` and `found` are, in the merge below.
        # Reusing them made every thread hammer one shared Set, which throws
        # "Multiple concurrent writes to Dict" and, before it throws, silently
        # cross-rejects other threads' signatures as duplicates.
        local_rng = Xoshiro(seeds[w])
        local_seen = Set{Vector{NTuple{4, Float64}}}()
        local_found = Floorplan[]
        for _ in 1:(40 * per_worker)
            (length(local_found) >= per_worker || time() > deadline) && break
            tries[w] += 1
            child = grow_nfs(pool[rand(local_rng, 1:length(pool))], local_rng;
                             max_repair=max_repair)
            child === nothing && continue
            sig = d4_signature(child)
            sig in local_seen && continue
            push!(local_seen, sig)
            push!(local_found, child)
        end
        parts[w] = local_found
    end
    merged = Floorplan[]
    seen = Set{Vector{NTuple{4, Float64}}}()
    for part in parts, fp in part
        length(merged) >= target && break
        sig = d4_signature(fp)
        sig in seen && continue
        push!(seen, sig)
        push!(merged, Floorplan("nfs_k$(length(fp))_$(length(merged) + 1)", fp.rectangles))
    end
    return merged, sum(tries)
end

function main()
    args = ARGS
    geti(flag, dflt) = (i = findfirst(==(flag), args);
                        i === nothing ? dflt : parse(Int, args[i + 1]))
    per_k    = geti("--per-k", 3000)
    kmax     = geti("--kmax", 18)
    seconds  = Float64(geti("--seconds", 120))

    cache = joinpath("data", "cache", "floorplans.cache")
    isfile(cache) || (println(stderr, "no cache at $cache; run init-cache"); return 1)
    lib = load_library(cache)
    println("corpus: $(length(lib)) cached child topologies")

    # How much of the EXISTING corpus is already face-share-free, by k.
    by_k = Dict{Int, Vector{Floorplan}}()
    for fp in lib
        push!(get!(by_k, length(fp), Floorplan[]), fp)
    end
    println("\n  k   cached   face-share-free   share")
    nfs_by_k = Dict{Int, Vector{Floorplan}}()
    for k in sort(collect(keys(by_k)))
        free = filter(is_face_share_free, by_k[k])
        nfs_by_k[k] = free
        println(rpad(k, 4), rpad(length(by_k[k]), 9),
                rpad(length(free), 18),
                string(round(100 * length(free) / length(by_k[k]), digits=1), "%"))
    end
    total_free = sum(length(v) for v in values(nfs_by_k))
    println("total face-share-free in corpus: $total_free / $(length(lib))")

    # Signatures of everything already cached, so "new" means new.
    println("\nsignature index of the existing corpus ...")
    known = Set{Vector{NTuple{4, Float64}}}()
    for fp in lib
        push!(known, d4_signature(fp))
    end

    println("\ngenerating, $(Threads.nthreads()) threads, $(seconds)s per k")
    println("  k   pool   attempts   found   new     rate")
    generated = Floorplan[]
    for k in 9:kmax
        pool = get(nfs_by_k, k - 1, Floorplan[])
        isempty(pool) && continue
        length(pool) > 20000 && (pool = pool[1:20000])
        level, tries = generate_level(pool, per_k; max_seconds=seconds)
        new = [fp for fp in level if !(d4_signature(fp) in known)]
        for fp in new
            push!(known, d4_signature(fp))
        end
        append!(generated, new)
        println(rpad(k, 4), rpad(length(pool), 7), rpad(tries, 11),
                rpad(length(level), 8), rpad(length(new), 8),
                string(round(100 * length(level) / max(tries, 1), digits=1), "%"))
    end

    out = joinpath("data", "import", "nfs_generated.cache")
    mkpath(dirname(out))
    save_library(out, generated)
    println("\nwrote $(length(generated)) new face-share-free topologies -> $out")
    return 0
end

main()
