#!/usr/bin/env julia
# Expand the OUTER ROOT pool -- the templates that split the square into R
# big boxes at the top of the hierarchy.
#
# The shipped pool has 339 templates for R = 2..14, and they are consumed in
# raw file order: solve() takes avail_roots[1:max_roots] with max_roots = 5
# for R = 7..10. So the entire top-level split of every large n is chosen
# from the first five entries of a list nobody ranked.
#
# The pool is also tiny in absolute terms. The number of rectangulations
# with R rectangles follows the Baxter numbers -- 24, 116, 642, 3938, 26194,
# 186042 for R = 4..9 -- so after dividing by the 8 symmetries of the square
# there are roughly 23,000 distinct R=9 templates. The cache holds 7.
#
#   julia -t auto tools/gen_roots.jl [--per-r N] [--rmax R] [--seconds S]
#
# Writes data/cache/outer_roots_expanded.cache, leaving the original alone.

include(joinpath(@__DIR__, "..", "src", "MondrianSolver.jl"))
using .MondrianSolver
using Random

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

"""One template with R+1 boxes from one with R: split, equalize, maybe flip.

No structural predicate is imposed. The point is coverage -- the ranking
inside solve() decides what is worth using, and it cannot rank what is not
there.
"""
function grow_root(fp::Floorplan, rng::AbstractRNG; flip_prob::Float64=0.5,
                   repair::Int=0)
    cand = split_rectangle(fp, rand(rng, 1:length(fp)), rand(rng, Bool),
                           0.2 + 0.6 * rand(rng))
    cur = equalize(cand)
    cur === nothing && return nothing

    # Splitting a rectangle ALWAYS leaves the two halves flush along their
    # shared edge, so a fresh split is never face-share-free. With repair=0
    # the walk is undirected and the pool comes out guillotine-flavoured:
    # of 48,838 templates generated that way, four were face-share-free
    # against 167 of the 339 shipped. With repair>0 the flips are retried
    # until the predicate is met, the same way grow_prime repairs toward
    # primality.
    if repair > 0
        for _ in 1:repair
            has_face_share(cur) || return cur
            flips = t_flips(cur)
            isempty(flips) && return nothing
            nxt = equalize(flips[rand(rng, 1:length(flips))])
            nxt === nothing && return nothing
            cur = nxt
        end
        return has_face_share(cur) ? nothing : cur
    end

    if rand(rng) < flip_prob
        flips = t_flips(cur)
        if !isempty(flips)
            nxt = equalize(flips[rand(rng, 1:length(flips))])
            nxt !== nothing && (cur = nxt)
        end
    end
    return cur
end

function grow_level(pool::Vector{Floorplan}, target::Int, R::Int;
                    max_seconds::Float64=60.0, repair::Int=0)
    isempty(pool) && return Floorplan[], 0
    nth = max(1, Threads.nthreads())
    deadline = time() + max_seconds
    seeds = [rand(UInt64) for _ in 1:nth]
    per_worker = cld(11 * target, 10 * nth)
    parts = [Floorplan[] for _ in 1:nth]
    tries = zeros(Int, nth)
    Threads.@threads for w in 1:nth
        local_rng = Xoshiro(seeds[w])
        local_seen = Set{Vector{NTuple{4, Float64}}}()
        local_found = Floorplan[]
        while length(local_found) < per_worker && time() < deadline
            tries[w] += 1
            child = grow_root(pool[rand(local_rng, 1:length(pool))], local_rng; repair=repair)
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
        push!(merged, Floorplan("root_R$(R)_$(length(merged) + 1)", fp.rectangles))
    end
    return merged, sum(tries)
end

function main()
    geti(f, d) = (i = findfirst(==(f), ARGS); i === nothing ? d : parse(Int, ARGS[i + 1]))
    per_r   = geti("--per-r", 4000)
    repair  = geti("--repair", 0)
    rmax    = geti("--rmax", 14)
    seconds = Float64(geti("--seconds", 60))

    base = repair > 0 ? "outer_roots_expanded.cache" : "outer_roots.cache"
    orig = load_library(joinpath("data", "cache", base))
    by_r = Dict{Int, Vector{Floorplan}}()
    for fp in orig
        push!(get!(by_r, length(fp), Floorplan[]), fp)
    end
    println("shipped pool: $(length(orig)) templates, R = "
            * string(minimum(keys(by_r))) * ".." * string(maximum(keys(by_r))))

    seen = Set{Vector{NTuple{4, Float64}}}()
    for fp in orig
        push!(seen, d4_signature(fp))
    end

    out = copy(orig)
    println("\n R   seed   attempts   kept    new     fsfree   total now")
    for R in 3:rmax
        pool = get(by_r, R - 1, Floorplan[])
        isempty(pool) && continue
        length(pool) > 8000 && (pool = pool[1:8000])
        level, tries = grow_level(pool, per_r, R; max_seconds=seconds, repair=repair)
        new = Floorplan[]
        for fp in level
            sig = d4_signature(fp)
            sig in seen && continue
            push!(seen, sig)
            push!(new, fp)
        end
        # Feed everything at this R forward, generated plus shipped, so the
        # next level is grown from the widest base available.
        by_r[R] = vcat(get(by_r, R, Floorplan[]), new)
        append!(out, new)
        fsf = count(fp -> !has_face_share(fp), new)
        println(rpad(R, 4), rpad(length(pool), 7), rpad(tries, 11), rpad(length(level), 8),
                rpad(length(new), 8), rpad(fsf, 9), length(by_r[R]))
    end

    dest = joinpath("data", "cache", "outer_roots_expanded.cache")
    println("seeding from $base, repair=$repair")
    save_library(dest, out)
    println("\nwrote $(length(out)) templates -> $dest "
            * "(was $(length(orig)))")
    return 0
end

main()
