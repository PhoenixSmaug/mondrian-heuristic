# ==============================================================================
# core.jl — Fundamental Data Models, Sentinels, & Geometric Validators
# ==============================================================================

# ─── 1. Constants & Primitive Aliases ────────────────────────────────────────

const LOW  = -1   # Sentinel for boundary x = 0 or y = 0
const HIGH = -2   # Sentinel for boundary x = 1 (or n) or y = 1 (or n)

const Rect    = NTuple{4, Float64}  # Normalized rectangle (x0, y0, x1, y1) in [0, 1]
const IntRect = NTuple{4, Int}      # Integer grid rectangle (x0, y0, x1, y1) in [0, n]
const RectRef = NTuple{4, Int}      # Rectangle walls (xl, y0, xr, y1) as wall IDs or sentinels
const Dim     = NTuple{2, Int}      # Canonical dimension (min(w, h), max(w, h))

@inline canonical_dim(w::Int, h::Int)::Dim = w <= h ? (w, h) : (h, w)

# ─── 2. Data Structures ──────────────────────────────────────────────────────

"""
    Floorplan(id, rectangles)

Geometric layout represented as non-overlapping rectangles tiling `[0, 1]²`.
"""
struct Floorplan
    id::String
    rectangles::Vector{Rect}
end

Base.length(fp::Floorplan) = length(fp.rectangles)

"""
    Topology(id, rects, position, is_x)

Combinatorial rectangulation defined by `v = k - 1` maximal internal wall segments.
- `rects`: 4-tuples `(xl, y0, xr, y1)` referencing wall indices or sentinels.
- `position`: Relative coordinates `p ∈ (0, 1)^v` of the walls.
- `is_x`: Boolean flags indicating vertical (`true`) or horizontal (`false`) walls.
"""
struct Topology
    id::String
    rects::Vector{RectRef}
    position::Vector{Float64}
    is_x::Vector{Bool}
end

Base.length(top::Topology) = length(top.rects)

"""
    Solution(n, defect, rectangles, method)

A validated integer partition of an `n × n` square into non-congruent rectangles.
"""
struct Solution
    n::Int
    defect::Int
    rectangles::Vector{IntRect}
    method::String
end

# ─── 3. Floorplan ↔ Topology Bridge ──────────────────────────────────────────

"""
    topology(fp::Floorplan) -> Topology

Extract maximal wall segment incidences from a normalized `Floorplan` using
a disjoint-set merge of collinear overlapping boundary segments.
"""
function topology(fp::Floorplan)::Topology
    rects = fp.rectangles
    k = length(rects)
    refs = [fill(0, 4) for _ in 1:k]
    initial = Float64[]
    orientation = Bool[]

    for (is_x, sides) in ((true, (1, 3)), (false, (2, 4)))
        by_pos = Dict{Float64, Vector{Tuple{Float64, Float64, Int, Int}}}()
        for (i, r) in enumerate(rects)
            lo, hi = is_x ? (r[2], r[4]) : (r[1], r[3])
            for side in sides
                coord = round(r[side], digits=9)
                push!(get!(by_pos, coord, Tuple{Float64,Float64,Int,Int}[]), (lo, hi, i, side))
            end
        end

        for (pos, occurrences) in by_pos
            if isapprox(pos, 0.0; atol=1e-7) || isapprox(pos, 1.0; atol=1e-7)
                sentinel = isapprox(pos, 0.0; atol=1e-7) ? LOW : HIGH
                for (_, _, rect_idx, side) in occurrences
                    refs[rect_idx][side] = sentinel
                end
                continue
            end

            # Merge overlapping collinear segments into maximal walls
            m = length(occurrences)
            parent = collect(1:m)
            find_root(a) = parent[a] == a ? a : (parent[a] = find_root(parent[a]))
            for a in 1:m, b in 1:a-1
                oa, ob = occurrences[a], occurrences[b]
                if max(oa[1], ob[1]) < min(oa[2], ob[2]) - 1e-9
                    ra, rb = find_root(a), find_root(b)
                    ra != rb && (parent[rb] = ra)
                end
            end

            groups = Dict{Int, Vector{Tuple{Float64, Float64, Int, Int}}}()
            for (j, occ) in enumerate(occurrences)
                push!(get!(groups, find_root(j), Tuple{Float64,Float64,Int,Int}[]), occ)
            end

            for group in values(groups)
                var_id = length(initial) + 1
                push!(initial, pos)
                push!(orientation, is_x)
                for (_, _, rect_idx, side) in group
                    refs[rect_idx][side] = var_id
                end
            end
        end
    end

    return Topology(fp.id, [Tuple(q) for q in refs], initial, orientation)
end

"""
    floorplan(top::Topology) -> Floorplan

Evaluate wall positions into normalized coordinates on `[0, 1]²`.
"""
function floorplan(top::Topology)::Floorplan
    p = top.position
    coord(idx) = idx == LOW ? 0.0 : idx == HIGH ? 1.0 : p[idx]
    rects = [(coord(r[1]), coord(r[2]), coord(r[3]), coord(r[4])) for r in top.rects]
    sort!(rects)
    return Floorplan(top.id, rects)
end

# ─── 4. Independent Validators ───────────────────────────────────────────────

"""
    is_valid(fp::Floorplan; tol=1e-6) -> Bool

Verify that a normalized floorplan forms a gapless, non-overlapping tiling of `[0, 1]²`.
"""
function is_valid(fp::Floorplan; tol::Float64=1e-6)::Bool
    total_area = 0.0
    rects = fp.rectangles
    for (i, r) in enumerate(rects)
        (r[1] < r[3] && r[2] < r[4]) || return false
        (r[1] >= -tol && r[2] >= -tol && r[3] <= 1.0 + tol && r[4] <= 1.0 + tol) || return false
        total_area += (r[3] - r[1]) * (r[4] - r[2])
        for j in 1:i-1
            s = rects[j]
            if max(r[1], s[1]) < min(r[3], s[3]) - tol && max(r[2], s[2]) < min(r[4], s[4]) - tol
                return false
            end
        end
    end
    return abs(total_area - 1.0) <= tol
end

"""
    validate_solution(n::Int, rects::Vector{IntRect}) -> Union{Int, Nothing}

Rigorously validate an integer Mondrian partition:
1. Integer bounds within `[0, n] × [0, n]`.
2. Exact area sum: `∑ Area_i = n²`.
3. Disjoint interiors (no overlaps).
4. Pairwise non-congruence: all canonical dimensions `(min(w,h), max(w,h))` are unique.

Returns `defect = max_area - min_area` if valid, or `nothing` if invalid.
"""
function validate_solution(n::Int, rects::Vector{IntRect})::Union{Int, Nothing}
    length(rects) >= 2 || return nothing
    dimensions = Set{Dim}()
    areas = Int[]

    for (i, r) in enumerate(rects)
        (0 <= r[1] < r[3] <= n && 0 <= r[2] < r[4] <= n) || return nothing
        w, h = r[3] - r[1], r[4] - r[2]
        d = canonical_dim(w, h)
        d in dimensions && return nothing
        push!(dimensions, d)
        push!(areas, w * h)

        for j in 1:i-1
            s = rects[j]
            if max(r[1], s[1]) < min(r[3], s[3]) && max(r[2], s[2]) < min(r[4], s[4])
                return nothing
            end
        end
    end

    sum(areas) == n * n || return nothing
    return maximum(areas) - minimum(areas)
end