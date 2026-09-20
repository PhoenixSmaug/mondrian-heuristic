# ==============================================================================
# geometry.jl — Continuous Newton Equalizer & Atomic Geometric Mutations
# ==============================================================================

# ─── 1. Continuous Newton-Raphson Area Continuation ──────────────────────────

@inline function wall_coord(idx::Int, p::Vector{Float64})::Float64
    idx == LOW ? 0.0 : idx == HIGH ? 1.0 : p[idx]
end

"""
    area_jacobian(top::Topology, p::Vector{Float64}) -> (areas, J)

Compute rectangle areas `A_i(p)` and the `(k-1) × (k-1)` analytic Jacobian matrix:
`J_{ij} = ∂A_i / ∂p_j`.
"""
function area_jacobian(top::Topology, p::Vector{Float64})
    k = length(top.rects)
    v = length(p)
    areas = zeros(Float64, k)
    J = zeros(Float64, k - 1, v)

    for (i, (xl, y0, xr, y1)) in enumerate(top.rects)
        w = wall_coord(xr, p) - wall_coord(xl, p)
        h = wall_coord(y1, p) - wall_coord(y0, p)
        areas[i] = w * h

        if i < k
            xl > 0 && (J[i, xl] -= h)
            xr > 0 && (J[i, xr] += h)
            y0 > 0 && (J[i, y0] -= w)
            y1 > 0 && (J[i, y1] += w)
        end
    end
    return areas, J
end

function is_positive_realization(top::Topology, p::Vector{Float64})::Bool
    all(0.0 < x < 1.0 for x in p) || return false
    for (xl, y0, xr, y1) in top.rects
        (wall_coord(xr, p) > wall_coord(xl, p) + 1e-9 &&
         wall_coord(y1, p) > wall_coord(y0, p) + 1e-9) || return false
    end
    return true
end

"""
    equalize_weighted(top::Topology, weights::Vector{Int}; tol=1e-12, max_iters=60) -> Union{Topology, Nothing}

Continue `top` to an exact geometric realization where each rectangle `i` has relative area
`weights[i] / ∑ weights`. Returns `nothing` if the topology is non-realizable or singular.
"""
function equalize_weighted(top::Topology, weights::Vector{Int};
                           tol::Float64=1e-12, max_iters::Int=60)::Union{Topology, Nothing}
    k = length(top.rects)
    length(weights) == k || error("weight count mismatch: expected $k, got $(length(weights))")
    length(top.position) == k - 1 || return nothing

    targets = Float64.(weights) ./ sum(weights)
    target_sub = targets[1:k-1]
    p = copy(top.position)

    for _ in 1:max_iters
        areas, J = area_jacobian(top, p)
        residual = areas[1:k-1] .- target_sub

        maximum(abs, residual) < tol && return Topology(top.id, top.rects, p, top.is_x)

        delta = try
            J \ (-residual)
        catch
            return nothing
        end

        old_norm = dot(residual, residual)
        alpha = 1.0
        accepted = false

        while alpha > 1e-10
            trial_p = p .+ alpha .* delta
            if is_positive_realization(top, trial_p)
                next_areas, _ = area_jacobian(top, trial_p)
                next_residual = next_areas[1:k-1] .- target_sub
                if dot(next_residual, next_residual) < old_norm
                    p = trial_p
                    accepted = true
                    break
                end
            end
            alpha *= 0.5
        end

        accepted || return nothing
    end
    return nothing
end

"""
    equalize(top::Topology; tol=1e-12, max_iters=60) -> Union{Topology, Nothing}

Continue `top` to an exact equal-area floorplan where every rectangle has area `1/k`.
"""
function equalize(top::Topology; tol::Float64=1e-12, max_iters::Int=60)::Union{Topology, Nothing}
    k = length(top.rects)
    return equalize_weighted(top, fill(1, k); tol=tol, max_iters=max_iters)
end

"""
    equalize(fp::Floorplan; tol=1e-12, max_iters=60) -> Union{Floorplan, Nothing}

Convenience overload to directly equalize a `Floorplan`.
"""
function equalize(fp::Floorplan; tol::Float64=1e-12, max_iters::Int=60)::Union{Floorplan, Nothing}
    top = topology(fp)
    eq_top = equalize(top; tol=tol, max_iters=max_iters)
    eq_top === nothing && return nothing
    return floorplan(eq_top)
end

# ─── 2. Atomic Geometric Mutation Operators ─────────────────────────────────

"""
    t_flips(fp::Floorplan) -> Vector{Floorplan}

Generate all elementary T-junction flips on adjacent rectangle pairs forming an L-shape.
Swaps the internal cut segment with the perpendicular reflex-corner cut.
"""
function t_flips(fp::Floorplan)::Vector{Floorplan}
    rects = fp.rectangles
    k = length(rects)
    candidates = Floorplan[]
    tol = 1e-6

    for i in 1:k, j in 1:i-1
        a, b = rects[i], rects[j]
        horiz = abs(a[4] - b[2]) < tol || abs(b[4] - a[2]) < tol
        vert  = abs(a[3] - b[1]) < tol || abs(b[3] - a[1]) < tol
        horiz || vert || continue

        replacement = nothing
        if horiz
            bot, top = a[4] <= b[2] + tol ? (a, b) : (b, a)
            ov_x0, ov_x1 = max(bot[1], top[1]), min(bot[3], top[3])
            if ov_x0 < ov_x1 - tol
                left_aligned  = abs(bot[1] - top[1]) < tol
                right_aligned = abs(bot[3] - top[3]) < tol
                if xor(left_aligned, right_aligned)
                    cut = left_aligned ? min(bot[3], top[3]) : max(bot[1], top[1])
                    wide = bot[3] > top[3] ? bot : top
                    if left_aligned
                        replacement = ((bot[1], bot[2], cut, top[4]),
                                       (cut, wide[2], max(bot[3], top[3]), wide[4]))
                    else
                        replacement = ((min(bot[1], top[1]), wide[2], cut, wide[4]),
                                       (cut, bot[2], bot[3], top[4]))
                    end
                end
            end
        elseif vert
            left, right = a[3] <= b[1] + tol ? (a, b) : (b, a)
            ov_y0, ov_y1 = max(left[2], right[2]), min(left[4], right[4])
            if ov_y0 < ov_y1 - tol
                bot_aligned = abs(left[2] - right[2]) < tol
                top_aligned = abs(left[4] - right[4]) < tol
                if xor(bot_aligned, top_aligned)
                    cut = bot_aligned ? min(left[4], right[4]) : max(left[2], right[2])
                    tall = left[4] > right[4] ? left : right
                    if bot_aligned
                        replacement = ((left[1], left[2], right[3], cut),
                                       (tall[1], cut, tall[3], max(left[4], right[4])))
                    else
                        replacement = ((tall[1], min(left[2], right[2]), tall[3], cut),
                                       (left[1], cut, right[3], left[4]))
                    end
                end
            end
        end

        if replacement !== nothing
            new_rects = [rects[z] for z in 1:k if z != i && z != j]
            push!(new_rects, replacement[1], replacement[2])
            sort!(new_rects)
            cand = Floorplan("$(fp.id)_flip", new_rects)
            is_valid(cand) && push!(candidates, cand)
        end
    end
    return candidates
end

function segment_spans(top::Topology)
    v = length(top.position)
    spans = fill((Inf, -Inf), v)
    for (i, isx) in enumerate(top.is_x)
        lo, hi = Inf, -Inf
        for q in top.rects
            if isx && (q[1] == i || q[3] == i)
                y0 = q[2] == LOW ? 0.0 : q[2] == HIGH ? 1.0 : top.position[q[2]]
                y1 = q[4] == LOW ? 0.0 : q[4] == HIGH ? 1.0 : top.position[q[4]]
                lo = min(lo, y0); hi = max(hi, y1)
            elseif !isx && (q[2] == i || q[4] == i)
                x0 = q[1] == LOW ? 0.0 : q[1] == HIGH ? 1.0 : top.position[q[1]]
                x1 = q[3] == LOW ? 0.0 : q[3] == HIGH ? 1.0 : top.position[q[3]]
                lo = min(lo, x0); hi = max(hi, x1)
            end
        end
        spans[i] = (lo, hi)
    end
    return spans
end

"""
    wall_slides(top::Topology) -> Vector{Topology}

Slide adjacent, oppositely directed T-junctions along shared host wall segments.
Creates non-slicing pinwheel mutations unreachable by simple T-flips.
"""
function wall_slides(top::Topology)::Vector{Topology}
    spans = segment_spans(top)
    candidates = Topology[]
    v = length(top.position)

    for host in 1:v
        host_is_x = top.is_x[host]
        host_coord = top.position[host]
        host_lo, host_hi = spans[host]

        events = Tuple{Float64, Int, Int}[]
        for moving in 1:v
            top.is_x[moving] == host_is_x && continue
            m_lo, m_hi = spans[moving]
            dir = 0
            if isapprox(m_hi, host_coord; atol=1e-5)
                dir = -1
            elseif isapprox(m_lo, host_coord; atol=1e-5)
                dir = 1
            end
            dir == 0 && continue
            along = top.position[moving]
            (host_lo + 1e-5 < along < host_hi - 1e-5) || continue
            push!(events, (along, moving, dir))
        end

        sort!(events, by=x -> x[1])
        for i in 1:length(events)-1
            a, b = events[i], events[i+1]
            a[3] == b[3] && continue
            p = copy(top.position)
            p[a[2]], p[b[2]] = p[b[2]], p[a[2]]

            if is_positive_realization(top, p)
                trial_top = Topology("$(top.id)_slide", top.rects, p, top.is_x)
                eq = equalize(trial_top)
                eq !== nothing && push!(candidates, eq)
            end
        end
    end
    return candidates
end

@inline function dihedral_transform(r::Rect, symmetry::Int)::Rect
    x0, y0, x1, y1 = r
    if symmetry >= 4
        x0, x1 = 1.0 - x1, 1.0 - x0
        symmetry -= 4
    end
    if symmetry == 0
        return (x0, y0, x1, y1)
    elseif symmetry == 1
        return (1.0 - y1, x0, 1.0 - y0, x1)
    elseif symmetry == 2
        return (1.0 - x1, 1.0 - y1, 1.0 - x0, 1.0 - y0)
    else
        return (y0, 1.0 - x1, y1, 1.0 - x0)
    end
end

"""Find all sub-rectangular blocks composed of 2 to `max_pieces` rectangles."""
function rectangular_blocks(fp::Floorplan; max_pieces::Int=8)
    rects = fp.rectangles
    k = length(rects)
    xs = sort!(unique(vcat([r[1] for r in rects], [r[3] for r in rects])))
    ys = sort!(unique(vcat([r[2] for r in rects], [r[4] for r in rects])))
    blocks = Tuple{Vector{Int}, Rect}[]

    for i1 in 1:length(xs)-1, j1 in 1:length(ys)-1
        x0, y0 = xs[i1], ys[j1]
        for i2 in i1+1:length(xs), j2 in j1+1:length(ys)
            x1, y1 = xs[i2], ys[j2]
            ids = Int[]
            area_sum = 0.0
            for (idx, r) in enumerate(rects)
                if x0 <= r[1] + 1e-7 && r[3] <= x1 + 1e-7 && y0 <= r[2] + 1e-7 && r[4] <= y1 + 1e-7
                    push!(ids, idx)
                    area_sum += (r[3] - r[1]) * (r[4] - r[2])
                end
            end
            if 2 <= length(ids) <= max_pieces && length(ids) < k
                if abs(area_sum - (x1 - x0) * (y1 - y0)) < 1e-6
                    push!(blocks, (ids, (x0, y0, x1, y1)))
                end
            end
        end
    end
    return blocks
end

"""Replace a sub-rectangular block with a rotated replacement topology."""
function substitute_block(fp::Floorplan, block_ids::Vector{Int}, box::Rect,
                          replacement::Floorplan, symmetry::Int)::Union{Floorplan, Nothing}
    bx0, by0, bx1, by1 = box
    bw, bh = bx1 - bx0, by1 - by0
    removed = Set(block_ids)
    new_rects = [r for (i, r) in enumerate(fp.rectangles) if !(i in removed)]

    for r in replacement.rectangles
        tr = dihedral_transform(r, symmetry)
        push!(new_rects, (bx0 + tr[1] * bw, by0 + tr[2] * bh,
                          bx0 + tr[3] * bw, by0 + tr[4] * bh))
    end
    sort!(new_rects)
    cand = Floorplan("$(fp.id)_sub$(length(replacement))", new_rects)
    is_valid(cand) || return nothing
    return cand
end

function mutate(fp::Floorplan; rng::AbstractRNG=Random.default_rng())::Union{Floorplan, Nothing}
    flips = t_flips(fp)
    isempty(flips) && return nothing
    return rand(rng, flips)
end