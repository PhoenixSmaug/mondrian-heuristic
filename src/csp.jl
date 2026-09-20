# ==============================================================================
# csp.jl — Discrete Grid Rounding, Box Caching, & Forward-Checking CSP
# ==============================================================================

# ─── 1. Single-Level Lattice Rounding ────────────────────────────────────────

function instantiate_square(top::Topology, n::Int, mask::Int)::Union{Vector{IntRect}, Nothing}
    v = length(top.position)
    p = Vector{Int}(undef, v)
    for i in 1:v
        base = floor(Int, top.position[i] * n)
        p[i] = base + ((mask >> (i - 1)) & 1)
    end
    coord(idx) = idx == LOW ? 0 : idx == HIGH ? n : p[idx]

    rects = IntRect[]
    for (xl, y0, xr, y1) in top.rects
        x0, y_0, x1, y_1 = coord(xl), coord(y0), coord(xr), coord(y1)
        (x0 < x1 && y_0 < y_1) || return nothing
        push!(rects, (x0, y_0, x1, y_1))
    end
    return rects
end

function evaluate_mask(top::Topology, n::Int, mask::Int)
    rects = instantiate_square(top, n, mask)
    rects === nothing && return (typemax(Int), nothing, IntRect[])

    dims = Set{Dim}()
    min_a, max_a = typemax(Int), 0
    duplicate_count = 0

    for r in rects
        w, h = r[3] - r[1], r[4] - r[2]
        d = canonical_dim(w, h)
        d in dims ? (duplicate_count += 1) : push!(dims, d)
        a = w * h
        min_a = min(min_a, a)
        max_a = max(max_a, a)
    end

    defect = max_a - min_a
    score = defect + duplicate_count * max(n, defect + 1)
    return (score, duplicate_count == 0 ? defect : nothing, rects)
end

function screen_topology(top::Topology, n::Int; max_samples::Int=256)
    v = length(top.position)
    v == 0 && return (typemax(Int), 0, IntRect[])

    nearest = 0
    for i in 1:v
        frac = top.position[i] * n - floor(top.position[i] * n)
        frac >= 0.5 && (nearest |= 1 << (i - 1))
    end

    best_score, best_defect, best_mask = typemax(Int), nothing, nearest
    best_rects = IntRect[]

    function check_mask(m::Int)
        score, defect, rects = evaluate_mask(top, n, m)
        if score < best_score
            best_score, best_defect, best_mask = score, defect, m
            best_rects = rects
        end
    end

    check_mask(nearest); check_mask(0); check_mask((1 << v) - 1)

    curr_mask = nearest
    for _ in 1:3
        improved = false
        curr_score, _, _ = evaluate_mask(top, n, curr_mask)
        for i in 0:v-1
            trial = xor(curr_mask, 1 << i)
            s, _, _ = evaluate_mask(top, n, trial)
            if s < curr_score
                curr_mask = trial
                curr_score = s
                improved = true
                check_mask(curr_mask)
            end
        end
        improved || break
    end

    modulus = 1 << min(v, 62)
    state = UInt64(hash((top.id, n), UInt(0x4d6f6e64)))
    for _ in 1:min(max_samples, modulus)
        state = state * UInt64(6364136223846793005) + UInt64(1442695040888963407)
        check_mask(Int(state % UInt64(modulus)))
    end

    return (best_defect, best_mask, best_rects)
end

# ─── 2. Local Candidate Evaluation & Box Caching ─────────────────────────────

struct LocalBoxCandidate
    score::Float64
    min_area::Int
    max_area::Int
    dims::Vector{Dim}
    local_rects::Vector{IntRect}
end

function fast_child_score(child::Topology, w::Int, h::Int, target::Float64, transposed::Bool)
    ew, eh = transposed ? (h, w) : (w, h)
    max_dev = 0.0
    dims = Dim[]

    for (xl, y0, xr, y1) in child.rects
        x0 = xl == LOW ? 0 : xl == HIGH ? ew : round(Int, child.position[xl] * (child.is_x[xl] ? ew : eh))
        x1 = xr == LOW ? 0 : xr == HIGH ? ew : round(Int, child.position[xr] * (child.is_x[xr] ? ew : eh))
        y0_c = y0 == LOW ? 0 : y0 == HIGH ? eh : round(Int, child.position[y0] * (child.is_x[y0] ? ew : eh))
        y1_c = y1 == LOW ? 0 : y1 == HIGH ? eh : round(Int, child.position[y1] * (child.is_x[y1] ? ew : eh))
        (x0 < x1 && y0_c < y1_c) || return Inf

        rw, rh = x1 - x0, y1_c - y0_c
        push!(dims, canonical_dim(rw, rh))
        max_dev = max(max_dev, abs(rw * rh - target))
    end

    length(unique(dims)) == length(dims) || return max_dev + (w * h)
    return max_dev
end

function enumerate_child_roundings(child::Topology, w::Int, h::Int, target::Float64, transposed::Bool;
                                   keep::Int=10)::Vector{LocalBoxCandidate}
    v = length(child.position)
    ew, eh = transposed ? (h, w) : (w, h)
    low_p = [floor(Int, child.position[i] * (child.is_x[i] ? ew : eh)) for i in 1:v]

    masks_to_test = if v <= 13
        0:(1 << v)-1
    else
        nearest = 0
        for i in 1:v
            (child.position[i] * (child.is_x[i] ? ew : eh) - low_p[i]) >= 0.5 && (nearest |= 1 << (i - 1))
        end
        m_set = Int[nearest, 0, (1 << v) - 1]
        for i in 0:v-1; push!(m_set, xor(nearest, 1 << i)); end
        rng_state = UInt64(0x12345678)
        for _ in 1:256
            rng_state = rng_state * UInt64(6364136223846793005) + UInt64(1442695040888963407)
            push!(m_set, Int(rng_state % UInt64(1 << min(v, 62))))
        end
        unique(m_set)
    end

    cands = LocalBoxCandidate[]
    p = Vector{Int}(undef, v)

    for mask in masks_to_test
        for i in 1:v; p[i] = low_p[i] + ((mask >> (i - 1)) & 1); end
        coord(idx, is_w) = idx == LOW ? 0 : idx == HIGH ? (is_w ? ew : eh) : p[idx]

        valid = true
        dims = Dim[]; rects = IntRect[]; areas = Int[]

        for (xl, y0, xr, y1) in child.rects
            rx0, ry0 = coord(xl, true), coord(y0, false)
            rx1, ry1 = coord(xr, true), coord(y1, false)
            (rx0 < rx1 && ry0 < ry1) || (valid = false; break)

            if transposed; rx0, ry0, rx1, ry1 = eh - ry1, rx0, eh - ry0, rx1; end
            rw, rh = rx1 - rx0, ry1 - ry0
            d = canonical_dim(rw, rh)
            d in dims && (valid = false; break)
            push!(dims, d); push!(areas, rw * rh)
            push!(rects, (rx0, ry0, rx1, ry1))
        end

        valid || continue
        min_a, max_a = minimum(areas), maximum(areas)
        score = max(abs(min_a - target), abs(max_a - target))
        push!(cands, LocalBoxCandidate(score, min_a, max_a, dims, rects))
    end

    sort!(cands, by=c -> c.score)
    return cands[1:min(length(cands), keep)]
end

function get_cached_box_candidates(children::Vector{Topology}, w::Int, h::Int, target::Float64,
                                   cache::Dict{Tuple{Int,Int}, Vector{LocalBoxCandidate}};
                                   shortlist_count::Int=64, keep_per_child::Int=8,
                                   cands_cap::Int=80)::Vector{LocalBoxCandidate}
    return get!(cache, (w, h)) do
        scored = Tuple{Float64, Int, Bool}[]
        for (idx, child) in enumerate(children)
            for transposed in (false, true)
                s = fast_child_score(child, w, h, target, transposed)
                isfinite(s) && push!(scored, (s, idx, transposed))
            end
        end
        sort!(scored, by=x -> x[1])

        all_cands = LocalBoxCandidate[]
        for (_, idx, transposed) in scored[1:min(length(scored), shortlist_count)]
            append!(all_cands, enumerate_child_roundings(children[idx], w, h, target, transposed; keep=keep_per_child))
        end

        sort!(all_cands, by=c -> c.score)
        unique_cands = LocalBoxCandidate[]
        seen = Set{Tuple{Int, Int, Vector{Dim}}}()
        for c in all_cands
            key = (c.min_area, c.max_area, sort(c.dims))
            if !(key in seen)
                push!(seen, key)
                push!(unique_cands, c)
                length(unique_cands) >= cands_cap && break
            end
        end
        unique_cands
    end
end

# ─── 3. Outer Mask Optimization ──────────────────────────────────────────────

function balanced_outer_masks(root::Topology, n::Int; keep::Int=20)
    v = length(root.position)
    candidates = Tuple{Int, Int, Vector{IntRect}}[]

    masks_to_check = if v <= 10
        0:(1 << v)-1
    else
        nearest = 0
        for i in 1:v
            (root.position[i] * n - floor(root.position[i] * n)) >= 0.5 && (nearest |= 1 << (i - 1))
        end
        m_set = Int[nearest, 0, (1 << v) - 1]
        for i in 0:v-1; push!(m_set, xor(nearest, 1 << i)); end
        rng_state = UInt64(0xabcdef01)
        for _ in 1:128
            rng_state = rng_state * UInt64(6364136223846793005) + UInt64(1442695040888963407)
            push!(m_set, Int(rng_state % UInt64(1 << min(v, 62))))
        end
        unique(m_set)
    end

    for mask in masks_to_check
        boxes = instantiate_square(root, n, mask)
        boxes === nothing && continue
        areas = [(b[3] - b[1]) * (b[4] - b[2]) for b in boxes]
        defect = maximum(areas) - minimum(areas)
        push!(candidates, (defect, mask, boxes))
    end

    sort!(candidates, by=x -> x[1])
    return candidates[1:min(length(candidates), keep)]
end

# ─── 4. Forward-Checking Interval CSP Engine ─────────────────────────────────

function solve_interval_csp(groups::Vector{Vector{LocalBoxCandidate}}, parent_boxes::Vector{IntRect};
                            best_defect_cap::Int=typemax(Int), max_proposals::Int=2500,
                            time_limit_sec::Float64=8.0)::Union{Vector{IntRect}, Nothing}
    any(isempty, groups) && return nothing
    R = length(groups)
    start_t = time()

    all_lows = sort!(unique(c.min_area for g in groups for c in g))
    proposals = Tuple{Int, Int, Int}[]

    for low in all_lows
        high = 0
        possible = true
        for g in groups
            local_high = typemax(Int)
            for c in g
                c.min_area >= low && (local_high = min(local_high, c.max_area))
            end
            if local_high == typemax(Int); possible = false; break; end
            high = max(high, local_high)
        end
        if possible
            base_span = high - low
            for slack in (0, 6, 12, 18, 24, 30)
                span = base_span + slack
                if span < best_defect_cap
                    push!(proposals, (span, low, high + slack))
                end
            end
        end
    end

    unique_proposals = unique(proposals)
    sort!(unique_proposals)
    max_search_props = min(length(unique_proposals), max_proposals)
    node_cap = R >= 8 ? 60_000 : 150_000

    for (_, low, high) in unique_proposals[1:max_search_props]
        (time() - start_t > time_limit_sec) && break

        eligible = [[c for c in g if c.min_area >= low && c.max_area <= high] for g in groups]
        any(isempty, eligible) && continue

        order = sortperm(1:R, by=i -> length(eligible[i]))
        chosen = Vector{LocalBoxCandidate}(undef, R)
        used_dims = Set{Dim}()
        nodes = 0

        function has_empty_future(depth::Int)::Bool
            for i in depth+1:R
                grp = order[i]
                viable = false
                for cand in eligible[grp]
                    if all(!(d in used_dims) for d in cand.dims)
                        viable = true
                        break
                    end
                end
                !viable && return true
            end
            return false
        end

        function dfs(depth::Int)::Bool
            nodes += 1
            if nodes > node_cap || (nodes % 256 == 0 && (time() - start_t > time_limit_sec))
                return false
            end
            depth > R && return true
            grp = order[depth]

            for cand in eligible[grp]
                if all(!(d in used_dims) for d in cand.dims)
                    for d in cand.dims; push!(used_dims, d); end
                    chosen[grp] = cand

                    # Forward checking lookahead
                    if depth < 2 || !has_empty_future(depth)
                        if dfs(depth + 1); return true; end
                    end

                    for d in cand.dims; delete!(used_dims, d); end
                end
            end
            return false
        end

        if dfs(1)
            all_rects = IntRect[]
            for (bx, cand) in zip(parent_boxes, chosen)
                bx0, by0 = bx[1], bx[2]
                for r in cand.local_rects
                    push!(all_rects, (bx0 + r[1], by0 + r[2], bx0 + r[3], by0 + r[4]))
                end
            end
            return all_rects
        end
    end
    return nothing
end