# ==============================================================================
# generators.jl — Planar Electrical Duals, Slicing Trees, & Cache Expansion
# ==============================================================================

# ─── 1. Planar Electrical Dual Network Generator ─────────────────────────────

function random_bipolar_faces(nv::Int, rng::AbstractRNG)
    nv >= 4 || error("triangulation requires at least 4 vertices")
    faces = Vector{Int}[ [1, 3, 2], [1, 2, 3] ]
    y_init = zeros(Float64, nv)
    y_init[1], y_init[2] = 0.0, 1.0
    y_init[3] = 0.3 + 0.4 * rand(rng)

    for v in 4:nv
        fi = rand(rng, 2:length(faces))
        a, b, c = faces[fi]
        lo = min(y_init[a], y_init[b], y_init[c])
        hi = max(y_init[a], y_init[b], y_init[c])
        y_init[v] = lo + (0.2 + 0.6 * rand(rng)) * (hi - lo)
        faces[fi] = [a, b, v]
        push!(faces, [b, c, v], [c, a, v])
    end
    return faces, y_init
end

function merge_faces(f1::Vector{Int}, f2::Vector{Int}, a::Int, b::Int)
    i1 = findfirst(i -> f1[i] == a && f1[mod1(i + 1, length(f1))] == b, 1:length(f1))
    i2 = findfirst(i -> f2[i] == b && f2[mod1(i + 1, length(f2))] == a, 1:length(f2))
    (i1 === nothing || i2 === nothing) && return nothing
    p1 = [f1[mod1(i1 + 1 + j, length(f1))] for j in 0:length(f1)-2]
    p2 = [f2[mod1(i2 + 1 + j, length(f2))] for j in 0:length(f2)-2]
    merged = vcat(p1, p2)
    length(unique(merged)) == length(merged) ? merged : nothing
end

function sparsify_faces!(faces::Vector{Vector{Int}}, y_init::Vector{Float64},
                         target_deletions::Int, rng::AbstractRNG)::Int
    deletions = 0
    for _ in 1:(40 * target_deletions)
        deletions >= target_deletions && break
        edge_map = Dict{Tuple{Int,Int},Int}()
        for (fi, f) in enumerate(faces)
            m = length(f)
            for j in 1:m
                edge_map[(f[j], f[mod1(j + 1, m)])] = fi
            end
        end

        candidates = [e for (e, f1) in edge_map if e[1] < e[2] && e != (1, 2) &&
                      f1 > 1 && get(edge_map, (e[2], e[1]), 1) > 1]
        isempty(candidates) && break
        a, b = rand(rng, candidates)
        f1, f2 = edge_map[(a, b)], edge_map[(b, a)]
        merged = merge_faces(faces[f1], faces[f2], a, b)
        merged === nothing && continue

        other_a = [v for f in faces for v in f if a in f && v != a && v != b]
        other_b = [v for f in faces for v in f if b in f && v != b && v != a]
        (a >= 3 && (!any(v -> y_init[v] < y_init[a], other_a) || !any(v -> y_init[v] > y_init[a], other_a))) && continue
        (b >= 3 && (!any(v -> y_init[v] < y_init[b], other_b) || !any(v -> y_init[v] > y_init[b], other_b))) && continue

        keep, drop = minmax(f1, f2)
        faces[keep] = merged
        deleteat!(faces, drop)
        deletions += 1
    end
    return deletions
end

function solve_electrical_potentials(nv::Int, sides::Dict{Tuple{Int,Int},Tuple{Int,Int}}, y_init::Vector{Float64})
    y = copy(y_init)
    internal = 3:nv
    m = length(internal)

    for _ in 1:60
        f = zeros(Float64, m)
        J = zeros(Float64, m, m)
        for (a, b) in keys(sides)
            u, v = y[a] < y[b] ? (a, b) : (b, a)
            gap = y[v] - y[u]
            inv_gap, inv_gap2 = 1.0 / gap, 1.0 / (gap * gap)
            if u >= 3
                ru = u - 2
                f[ru] -= inv_gap
                J[ru, ru] -= inv_gap2
                v >= 3 && (J[ru, v - 2] += inv_gap2)
            end
            if v >= 3
                rv = v - 2
                f[rv] += inv_gap
                J[rv, rv] -= inv_gap2
                u >= 3 && (J[rv, u - 2] += inv_gap2)
            end
        end

        maximum(abs, f) < 1e-10 && break
        step = J \ (-f)

        alpha = 1.0
        accepted = false
        for _ in 1:25
            trial_y = copy(y)
            for (i, v) in enumerate(internal)
                trial_y[v] += alpha * step[i]
            end
            if all(begin
                       u, v = y_init[a] < y_init[b] ? (a, b) : (b, a)
                       trial_y[v] > trial_y[u] + 1e-8
                   end for (a, b) in keys(sides)) &&
               all(0.0 < trial_y[v] < 1.0 for v in internal)
                y = trial_y
                accepted = true
                break
            end
            alpha *= 0.5
        end
        accepted || error("electrical Newton solve failed to converge")
    end
    return y
end

function solve_dual_potentials(nf::Int, sides::Dict{Tuple{Int,Int},Tuple{Int,Int}},
                               y::Vector{Float64}, st_left::Int)
    x = fill(NaN, nf)
    x[st_left] = 0.0
    adj = [Tuple{Int,Float64}[] for _ in 1:nf]

    for ((a, b), (minleft, minright)) in sides
        (a, b) == (1, 2) && continue
        lo, hi = y[a] < y[b] ? (a, b) : (b, a)
        left, right = lo == a ? (minleft, minright) : (minright, minleft)
        q = 1.0 / (y[hi] - y[lo])
        push!(adj[left], (right, q))
        push!(adj[right], (left, -q))
    end

    queue = Int[st_left]
    while !isempty(queue)
        f = popfirst!(queue)
        for (g, delta) in adj[f]
            if isnan(x[g])
                x[g] = x[f] + delta
                push!(queue, g)
            end
        end
    end

    all(isfinite, x) || error("dual graph disconnected")
    xmin, xmax = extrema(x)
    return (x .- xmin) ./ (xmax - xmin)
end

"""
    generate_electrical(k::Int; rng=default_rng(), id="") -> Floorplan

Construct an exact equal-area electrical dual floorplan with EXACTLY `k` pieces.
Uses random triangulations with edge sparsification (k = 3*nv - 7 - deletions).
"""
function generate_electrical(k::Int; rng::AbstractRNG=Random.default_rng(), id::String="")::Floorplan
    for _ in 1:25
        extra_v = rand(rng, 0:2)
        nv = max(4, ceil(Int, (k + 7) / 3) + extra_v)
        deletions = (3 * nv - 7) - k
        deletions < 0 && continue

        faces, y_init = random_bipolar_faces(nv, rng)
        if deletions > 0
            actual_del = sparsify_faces!(faces, y_init, deletions, rng)
            actual_del == deletions || continue
        end

        sides = Dict{Tuple{Int,Int},Tuple{Int,Int}}()
        for (fi, f) in enumerate(faces)
            m = length(f)
            for j in 1:m
                a, b = f[j], f[mod1(j + 1, m)]
                edge = minmax(a, b)
                cur = get(sides, edge, (0, 0))
                sides[edge] = a < b ? (fi, cur[2]) : (cur[1], fi)
            end
        end

        haskey(sides, (1, 2)) || continue

        y = try solve_electrical_potentials(nv, sides, y_init) catch; continue end
        st_left = sides[(1, 2)][1]
        x = try solve_dual_potentials(length(faces), sides, y, st_left) catch; continue end

        rectangles = Rect[]
        for ((a, b), (minleft, minright)) in sides
            (a, b) == (1, 2) && continue
            lo, hi = y[a] < y[b] ? (a, b) : (b, a)
            left, right = lo == a ? (minleft, minright) : (minright, minleft)
            x0, x1 = minmax(x[left], x[right])
            push!(rectangles, (x0, y[lo], x1, y[hi]))
        end

        length(rectangles) == k || continue
        sort!(rectangles)
        fp_id = isempty(id) ? "electrical_k$(k)" : id
        fp = Floorplan(fp_id, rectangles)
        is_valid(fp) && return fp
    end

    return generate_slicing(k; rng=rng, id=id)
end

function generate_electrical_library(count::Int; k_range::UnitRange{Int}=10:15,
                                     rng::AbstractRNG=Random.default_rng())::Vector{Floorplan}
    library = Floorplan[]
    seen = Set{UInt64}()
    attempts = 0
    max_attempts = 40 * count

    while length(library) < count && attempts < max_attempts
        attempts += 1
        k = rand(rng, k_range)
        fp = generate_electrical(k; rng=rng)
        h = hash(fp.rectangles)
        if !(h in seen) && is_valid(fp)
            push!(seen, h)
            push!(library, fp)
        end
    end
    return library
end

# ─── 2. Slicing Tree Generator ───────────────────────────────────────────────

function generate_slicing(k::Int; rng::AbstractRNG=Random.default_rng(), id::String="")::Floorplan
    k >= 2 || error("slicing floorplan requires at least 2 pieces")
    rects = Rect[]

    function descend(x0::Float64, y0::Float64, x1::Float64, y1::Float64, leaves::Int)
        if leaves == 1
            push!(rects, (x0, y0, x1, y1))
            return
        end
        a = rand(rng, 1:leaves-1)
        w, h = x1 - x0, y1 - y0

        cut_vertical = w > 1.15h ? rand(rng) < 0.85 :
                       h > 1.15w ? rand(rng) < 0.15 : rand(rng, Bool)

        if cut_vertical
            cut = x0 + w * (a / leaves)
            descend(x0, y0, cut, y1, a)
            descend(cut, y0, x1, y1, leaves - a)
        else
            cut = y0 + h * (a / leaves)
            descend(x0, y0, x1, cut, a)
            descend(x0, cut, x1, y1, leaves - a)
        end
    end

    descend(0.0, 0.0, 1.0, 1.0, k)
    sort!(rects)
    fp_id = isempty(id) ? "slicing_k$(k)" : id
    return Floorplan(fp_id, rects)
end

function generate_slicing_library(count::Int; k_range::UnitRange{Int}=10:16,
                                  rng::AbstractRNG=Random.default_rng())::Vector{Floorplan}
    library = Floorplan[]
    seen = Set{UInt64}()
    attempts = 0
    max_attempts = 150 * count

    while length(library) < count && attempts < max_attempts
        attempts += 1
        k = rand(rng, k_range)
        fp = generate_slicing(k; rng=rng)
        h = hash(fp.rectangles)
        if !(h in seen) && is_valid(fp)
            push!(seen, h)
            push!(library, fp)
        end
    end
    return library
end

# ─── 3. Full-Scale Cache Expansion Engines ───────────────────────────────────

"""
    expand_electrical_portfolio(base_by_k, target_per_k; rng) -> Vector{Floorplan}

Expand base electrical dual topologies up to `target_per_k` (e.g. 20,000)
per piece count k using T-flips and wall slides.
"""
function expand_electrical_portfolio(base_by_k::Dict{Int, Vector{Floorplan}}, target_per_k::Int;
                                     rng::AbstractRNG=Random.default_rng(),
                                     on_step=nothing)::Vector{Floorplan}
    expanded = Floorplan[]
    for k in sort(collect(keys(base_by_k)))
        base = base_by_k[k]
        target = max(length(base), target_per_k)
        start_t = time()

        archive = copy(base)
        seen = Set{Vector{Tuple{Float64,Float64,Float64,Float64}}}()
        for fp in archive
            push!(seen, [(round(r[1]; digits=6), round(r[2]; digits=6),
                          round(r[3]; digits=6), round(r[4]; digits=6)) for r in fp.rectangles])
        end

        frontier = copy(base)
        while length(archive) < target && !isempty(frontier)
            next_frontier = Floorplan[]
            for parent in frontier[randperm(rng, length(frontier))]
                length(archive) >= target && break

                # T-Flips
                for cand in t_flips(parent)
                    length(archive) >= target && break
                    eq = equalize(cand)
                    eq === nothing && continue
                    sig = [(round(r[1]; digits=6), round(r[2]; digits=6),
                            round(r[3]; digits=6), round(r[4]; digits=6)) for r in eq.rectangles]
                    if !(sig in seen)
                        push!(seen, sig); push!(archive, eq); push!(next_frontier, eq)
                    end
                end

                # Wall Slides
                top = topology(parent)
                for s in wall_slides(top)
                    length(archive) >= target && break
                    cand = floorplan(s)
                    sig = [(round(r[1]; digits=6), round(r[2]; digits=6),
                            round(r[3]; digits=6), round(r[4]; digits=6)) for r in cand.rectangles]
                    if !(sig in seen)
                        push!(seen, sig); push!(archive, cand); push!(next_frontier, cand)
                    end
                end
            end
            frontier = next_frontier
        end
        append!(expanded, archive)
        on_step !== nothing && on_step(k, length(base), length(archive), time() - start_t)
    end
    return expanded
end

"""
    generate_quilt_portfolio(small_fps, target_count; mink=10, maxk=18, rng) -> Vector{Floorplan}

Multi-level quilt expansion using dihedral block substitutions, T-flips, and
wall slides on `oeis_solutions.txt` seeds.
"""
function generate_quilt_portfolio(small_fps::Vector{Floorplan}, target_count::Int;
                                  mink::Int=10, maxk::Int=18,
                                  rng::AbstractRNG=Random.default_rng(),
                                  on_step=nothing)::Vector{Floorplan}
    libraries_by_m = Dict{Int, Vector{Floorplan}}()
    for fp in small_fps
        m = length(fp)
        2 <= m <= 8 && push!(get!(libraries_by_m, m, Floorplan[]), fp)
    end

    base = [fp for fp in small_fps if mink <= length(fp) <= maxk]
    isempty(base) && return Floorplan[]

    archive = copy(base)
    seen = Set{Vector{Tuple{Float64,Float64,Float64,Float64}}}()
    for fp in archive
        push!(seen, [(round(r[1]; digits=6), round(r[2]; digits=6),
                      round(r[3]; digits=6), round(r[4]; digits=6)) for r in fp.rectangles])
    end

    frontier = copy(base)
    depth = 0
    start_t = time()

    while length(archive) < target_count && !isempty(frontier)
        depth += 1
        next_frontier = Floorplan[]

        for parent in frontier[randperm(rng, length(frontier))]
            length(archive) >= target_count && break

            # 1. T-Flips
            for cand in t_flips(parent)
                length(archive) >= target_count && break
                eq = equalize(cand)
                eq === nothing && continue
                sig = [(round(r[1]; digits=6), round(r[2]; digits=6),
                        round(r[3]; digits=6), round(r[4]; digits=6)) for r in eq.rectangles]
                if !(sig in seen)
                    push!(seen, sig); push!(archive, eq); push!(next_frontier, eq)
                end
            end

            # 2. Wall Slides
            top = topology(parent)
            for s in wall_slides(top)
                length(archive) >= target_count && break
                cand = floorplan(s)
                sig = [(round(r[1]; digits=6), round(r[2]; digits=6),
                        round(r[3]; digits=6), round(r[4]; digits=6)) for r in cand.rectangles]
                if !(sig in seen)
                    push!(seen, sig); push!(archive, cand); push!(next_frontier, cand)
                end
            end

            # 3. Dihedral Block Substitutions
            blocks = rectangular_blocks(parent; max_pieces=6)
            for (ids, box) in blocks
                m = length(ids)
                haskey(libraries_by_m, m) || continue
                for repl in libraries_by_m[m][1:min(length(libraries_by_m[m]), 3)]
                    sym = rand(rng, 0:7)
                    cand = substitute_block(parent, ids, box, repl, sym)
                    cand === nothing && continue
                    eq = equalize(cand)
                    eq === nothing && continue
                    sig = [(round(r[1]; digits=6), round(r[2]; digits=6),
                            round(r[3]; digits=6), round(r[4]; digits=6)) for r in eq.rectangles]
                    if !(sig in seen)
                        push!(seen, sig); push!(archive, eq); push!(next_frontier, eq)
                    end
                    length(archive) >= target_count && break
                end
                length(archive) >= target_count && break
            end
        end

        on_step !== nothing && on_step(depth, length(archive), target_count, time() - start_t)
        frontier = next_frontier
    end
    return archive
end