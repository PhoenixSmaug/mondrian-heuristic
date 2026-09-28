#!/usr/bin/env julia
# Stock the OUTER ROOT pool with face-share-free templates from the prime
# catalogue already on disk.
#
# A prime rectangulation cannot contain a full-span tile -- the complement of
# such a tile would itself be a block of k-1 rectangles -- so every prime is
# face-share-free in the strong sense, boundary included. The repository
# already holds ~253,000 of them: the complete A340984 catalogue for k <= 14
# in data/import/a340984, and generated classes for k = 15..20 in
# data/import/prime_generated.cache. None of it has ever been used as outer
# root material, where the shipped pool is 339 templates total.
#
# These go in as OPTIONS. Nothing here prefers them; the ranking inside
# solve() decides on merit, and it cannot pick what is absent.
#
#   julia -t auto tools/roots_from_primes.jl [--max-per-r N] [--rmax R]

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

function main()
    geti(f, d) = (i = findfirst(==(f), ARGS); i === nothing ? d : parse(Int, ARGS[i + 1]))
    max_per_r = geti("--max-per-r", 12000)
    rmax      = geti("--rmax", 20)

    base = joinpath("data", "cache", "outer_roots_expanded.cache")
    isfile(base) || (base = joinpath("data", "cache", "outer_roots.cache"))
    pool = load_library(base)
    println("existing pool: $(length(pool)) from $base")

    seen = Set{Vector{NTuple{4, Float64}}}()
    for fp in pool
        push!(seen, d4_signature(fp))
    end

    cands = Floorplan[]
    cat = joinpath("data", "import", "a340984")
    if isdir(cat)
        c = parse_a340984(cat; k_range = 2:min(14, rmax), max_per_k = max_per_r)
        println("A340984 catalogue: $(length(c)) classes")
        append!(cands, c)
    end
    gen = joinpath("data", "import", "prime_generated.cache")
    if isfile(gen) && rmax > 14
        g = load_library(gen)
        g = filter(fp -> 15 <= length(fp) <= rmax, g)
        println("generated primes k=15..$rmax: $(length(g))")
        append!(cands, g)
    end
    isempty(cands) && (println(stderr, "no prime sources found"); return 1)

    # The catalogue stores integer drawings normalised onto the unit square.
    # Outer roots are instantiated from equalized wall positions, so each one
    # has to be run through the Newton equalizer; the ones that do not
    # converge to an equal-area realization cannot serve as roots.
    println("\nequalizing $(length(cands)) candidates on $(Threads.nthreads()) threads ...")
    eq = Vector{Union{Floorplan, Nothing}}(undef, length(cands))
    Threads.@threads for i in eachindex(cands)
        eq[i] = try
            equalize(cands[i])
        catch
            nothing
        end
    end

    kept  = Dict{Int, Int}()
    added = Floorplan[]
    nfs   = Dict{Int, Int}()
    for fp in eq
        fp === nothing && continue
        R = length(fp)
        (2 <= R <= rmax) || continue
        get(kept, R, 0) >= max_per_r && continue
        sig = d4_signature(fp)
        sig in seen && continue
        push!(seen, sig)
        kept[R] = get(kept, R, 0) + 1
        has_face_share(fp) || (nfs[R] = get(nfs, R, 0) + 1)
        push!(added, Floorplan("prime_root_R$(R)_$(get(kept, R, 0))", fp.rectangles))
    end

    by_r = Dict{Int, Int}()
    for fp in pool
        by_r[length(fp)] = get(by_r, length(fp), 0) + 1
    end
    println("\n R   had     added   of which fsfree   total")
    for R in 2:rmax
        (get(by_r, R, 0) == 0 && get(kept, R, 0) == 0) && continue
        println(rpad(R, 4), rpad(get(by_r, R, 0), 8), rpad(get(kept, R, 0), 8),
                rpad(get(nfs, R, 0), 18), get(by_r, R, 0) + get(kept, R, 0))
    end

    append!(pool, added)
    dest = joinpath("data", "cache", "outer_roots_expanded.cache")
    save_library(dest, pool)
    println("\nwrote $(length(pool)) templates -> $dest (added $(length(added)))")
    return 0
end

main()
