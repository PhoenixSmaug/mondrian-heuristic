# ==============================================================================
# prime.jl — Prime (Block-Irreducible) Rectangulations: Test & Generation
# ==============================================================================
#
# A rectangulation is *prime* (A340984) when no sub-collection of 2 .. k-1 of its
# rectangles has a rectangular union. In particular no two rectangles share a
# complete face, and the tiling cannot be assembled by substituting a smaller
# rectangulation into a rectangle of another one — it has to be found directly.
# OEIS A340984 catalogues every class up to k = 14; beyond that this file grows
# new ones from the catalogue by a split-and-repair random walk.

# ─── 1. Primality Test ───────────────────────────────────────────────────────

"""
    coarse_grid(fp) -> (boxes, owner, nx, ny)

Index a floorplan onto the grid of its distinct wall coordinates. `boxes[i]` is
rectangle `i` in grid indices, `owner[a, c]` is the rectangle covering cell
`(a, c)`.
"""
function coarse_grid(fp::Floorplan)
    rects = fp.rectangles
    xs = sort!(unique(round.(vcat([r[1] for r in rects], [r[3] for r in rects]); digits=9)))
    ys = sort!(unique(round.(vcat([r[2] for r in rects], [r[4] for r in rects]); digits=9)))
    ix(v) = searchsortedfirst(xs, round(v; digits=9))
    iy(v) = searchsortedfirst(ys, round(v; digits=9))

    boxes = [(ix(r[1]), iy(r[2]), ix(r[3]), iy(r[4])) for r in rects]
    owner = zeros(Int, length(xs) - 1, length(ys) - 1)
    for (i, b) in enumerate(boxes)
        for a in b[1]:b[3]-1, c in b[2]:b[4]-1
            owner[a, c] = i
        end
    end
    return boxes, owner, length(xs), length(ys)
end

"""
    has_rect_block(fp::Floorplan) -> Bool

True when some proper sub-collection of 2 .. k-1 rectangles tiles a rectangle.

A cheap `O(k²)` scan for full-face-sharing pairs runs first, which rejects the
overwhelming majority of candidates. The exhaustive pass then walks every grid
aligned box; a rectangle meeting a box without being contained in it must cross
the box boundary, so only the box perimeter has to be inspected.
"""
function has_rect_block(fp::Floorplan)::Bool
    rects = fp.rectangles
    k = length(rects)
    k >= 3 || return false

    for i in 1:k, j in 1:i-1
        a, b = rects[i], rects[j]
        (a[1] == b[1] && a[3] == b[3] && (a[4] == b[2] || b[4] == a[2])) && return true
        (a[2] == b[2] && a[4] == b[4] && (a[3] == b[1] || b[3] == a[1])) && return true
    end

    boxes, owner, nx, ny = coarse_grid(fp)
    @inbounds for i1 in 1:nx-1, j1 in 1:ny-1, i2 in i1+1:nx, j2 in j1+1:ny
        (i2 - i1 == nx - 1 && j2 - j1 == ny - 1) && continue   # the whole square

        closed = true
        for a in i1:i2-1, c in j1:j2-1
            (a == i1 || a == i2 - 1 || c == j1 || c == j2 - 1) || continue
            b = boxes[owner[a, c]]
            if !(b[1] >= i1 && b[3] <= i2 && b[2] >= j1 && b[4] <= j2)
                closed = false
                break
            end
        end
        closed || continue

        inside = 0
        for b in boxes
            (b[1] >= i1 && b[3] <= i2 && b[2] >= j1 && b[4] <= j2) && (inside += 1)
        end
        2 <= inside < k && return true
    end
    return false
end

"""
    is_prime_rectangulation(fp::Floorplan) -> Bool

True when `fp` is a prime rectangulation in the sense of OEIS A340984.
"""
is_prime_rectangulation(fp::Floorplan)::Bool = !has_rect_block(fp)

"""
    d4_signature(fp::Floorplan) -> Vector{NTuple{4, Float64}}

Orientation-independent key for an equal-area floorplan: the lexicographically
least rounded rectangle list over the eight dihedral images. The equal-area
realization of a generic rectangulation is unique, so this identifies the
topology up to symmetry exactly as A340984's `d4_certificate` does.
"""
function d4_signature(fp::Floorplan)::Vector{NTuple{4, Float64}}
    best = nothing
    for s in 0:7
        sig = sort!([round.(dihedral_transform(r, s); digits=6) for r in fp.rectangles])
        (best === nothing || sig < best) && (best = sig)
    end
    return best
end

"""
    grid_drawing(fp::Floorplan) -> (width, height, Vector{IntRect})

Snap a floorplan onto the integer grid of its distinct wall coordinates. This is
the drawing convention of the A340984 catalogue, so the result can be written
back out in the same schema and rendered by its `render.py`.
"""
function grid_drawing(fp::Floorplan)
    rects = fp.rectangles
    xs = sort!(unique(round.(vcat([r[1] for r in rects], [r[3] for r in rects]); digits=9)))
    ys = sort!(unique(round.(vcat([r[2] for r in rects], [r[4] for r in rects]); digits=9)))
    ix(v) = searchsortedfirst(xs, round(v; digits=9)) - 1
    iy(v) = searchsortedfirst(ys, round(v; digits=9)) - 1
    drawn = IntRect[(ix(r[1]), iy(r[2]), ix(r[3]), iy(r[4])) for r in rects]
    return length(xs) - 1, length(ys) - 1, drawn
end

# ─── 2. Split-and-Repair Growth ──────────────────────────────────────────────

"""
    split_rectangle(fp, i, vertical, frac) -> Floorplan

Cut rectangle `i` in two, giving a `k + 1` piece rectangulation. The result is
never prime — the two halves form a block — so it has to be repaired.
"""
function split_rectangle(fp::Floorplan, i::Int, vertical::Bool, frac::Float64)::Floorplan
    r = fp.rectangles[i]
    rects = [fp.rectangles[j] for j in eachindex(fp.rectangles) if j != i]
    if vertical
        c = r[1] + frac * (r[3] - r[1])
        push!(rects, (r[1], r[2], c, r[4]), (c, r[2], r[3], r[4]))
    else
        c = r[2] + frac * (r[4] - r[2])
        push!(rects, (r[1], r[2], r[3], c), (r[1], c, r[3], r[4]))
    end
    sort!(rects)
    return Floorplan(fp.id, rects)
end

"""
    grow_prime(fp, rng; max_repair=4, id="") -> Union{Floorplan, Nothing}

Split a random rectangle of the prime floorplan `fp` and try to repair the
resulting block away with up to `max_repair` random T-flips, re-equalizing after
every move. Returns an equal-area prime floorplan with one more rectangle, or
`nothing` if this attempt failed.

Short repair walks are deliberate: most successes come from one or two flips, so
a small `max_repair` maximizes throughput.
"""
function grow_prime(fp::Floorplan, rng::AbstractRNG;
                    max_repair::Int=4, id::String="")::Union{Floorplan, Nothing}
    cand = split_rectangle(fp, rand(rng, 1:length(fp)), rand(rng, Bool), 0.25 + 0.5 * rand(rng))
    cur = equalize(cand)
    cur === nothing && return nothing

    for _ in 1:max_repair
        if is_prime_rectangulation(cur)
            fid = isempty(id) ? "prime_gen_k$(length(cur))" : id
            return Floorplan(fid, cur.rectangles)
        end
        flips = t_flips(cur)
        isempty(flips) && return nothing
        nxt = equalize(flips[rand(rng, 1:length(flips))])
        nxt === nothing && return nothing
        cur = nxt
    end
    return nothing
end

# ─── 3. Level and Ladder Generators ──────────────────────────────────────────

"""
    generate_prime_level(pool, target; rng, max_seconds, max_repair, nthreads) -> Vector{Floorplan}

Grow `target` distinct prime rectangulations with `k + 1` pieces from a `pool` of
prime `k` piece parents. Duplicates are removed up to the dihedral group, so the
result is a set of distinct A340984 classes. Work is split across `nthreads`
workers, each with its own stream; the streams are merged and de-duplicated once
at the end.
"""
function generate_prime_level(pool::Vector{Floorplan}, target::Int;
                              rng::AbstractRNG=Random.default_rng(),
                              max_seconds::Float64=600.0,
                              max_repair::Int=4,
                              nthreads::Int=Threads.nthreads())::Vector{Floorplan}
    (isempty(pool) || target <= 0) && return Floorplan[]
    nthreads = max(1, min(nthreads, Threads.nthreads()))
    deadline = time() + max_seconds
    seeds = [rand(rng, UInt64) for _ in 1:nthreads]
    # 10% overshoot absorbs the duplicates found across independent workers.
    per_worker = cld(11 * target, 10 * nthreads)

    parts = [Floorplan[] for _ in 1:nthreads]
    Threads.@threads for w in 1:nthreads
        local_rng = Xoshiro(seeds[w])
        local_seen = Set{Vector{NTuple{4, Float64}}}()
        found = Floorplan[]
        # About 6% of attempts succeed, so 40 per wanted tiling is >2x headroom.
        budget = 40 * per_worker
        for _ in 1:budget
            length(found) >= per_worker && break
            (time() > deadline) && break
            child = grow_prime(pool[rand(local_rng, 1:length(pool))], local_rng; max_repair=max_repair)
            child === nothing && continue
            sig = d4_signature(child)
            sig in local_seen && continue
            push!(local_seen, sig)
            push!(found, child)
        end
        parts[w] = found
    end

    merged = Floorplan[]
    seen = Set{Vector{NTuple{4, Float64}}}()
    for part in parts, fp in part
        length(merged) >= target && break
        sig = d4_signature(fp)
        sig in seen && continue
        push!(seen, sig)
        push!(merged, Floorplan("prime_gen_k$(length(fp))_$(length(merged) + 1)", fp.rectangles))
    end
    return merged
end

"""
    generate_prime_ladder(base, kmax, target_per_k; rng, seconds_per_k, max_repair, on_step)

Climb from the largest available catalogue level up to `kmax`, generating
`target_per_k` distinct prime rectangulations at each new piece count and using
them as the parents for the next. Returns the generated levels only; `base` (the
complete A340984 catalogue) is left untouched.
"""
function generate_prime_ladder(base::Vector{Floorplan}, kmax::Int, target_per_k::Int;
                               rng::AbstractRNG=Random.default_rng(),
                               seconds_per_k::Float64=600.0,
                               max_repair::Int=4,
                               nthreads::Int=Threads.nthreads(),
                               on_step=nothing)::Dict{Int, Vector{Floorplan}}
    levels = Dict{Int, Vector{Floorplan}}()
    isempty(base) && return levels

    by_k = Dict{Int, Vector{Floorplan}}()
    for fp in base
        push!(get!(by_k, length(fp), Floorplan[]), fp)
    end
    k0 = maximum(keys(by_k))
    pool = by_k[k0]

    for k in (k0 + 1):kmax
        t0 = time()
        level = generate_prime_level(pool, target_per_k; rng=rng, max_seconds=seconds_per_k,
                                     max_repair=max_repair, nthreads=nthreads)
        on_step !== nothing && on_step(k, length(pool), length(level), time() - t0)
        isempty(level) && break
        levels[k] = level
        pool = level
    end
    return levels
end
