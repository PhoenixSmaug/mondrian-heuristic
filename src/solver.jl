# ==============================================================================
# solver.jl — Adaptive Strategy & Master Solve Coordination
# ==============================================================================
#
# The decision layer: every choice ABOUT the machinery in geometry.jl and
# csp.jl -- which hierarchy shapes (R outer boxes, each split into M child
# rectangles) to try, which outer roots and rounding masks to spend CSP calls
# on, and how to share a fixed budget of CSP calls between them.
#
# This layer was evolved with ShinkaEvolve, an LLM-driven program-evolution
# loop (proposals written by Claude), scoring each candidate by 50 x mean
# d_best_known / d_found over 30 values of n in 300..3000 and validating on 50
# disjoint held-out n. Two of the evolved programs are kept because they are
# architecturally different and win on different n; they share everything
# except the scheduler that allocates the CSP-call budget:
#
#   schedule = :bandit  successive halving with widening ROLLED INTO the
#                       halving loop (run_incremental_bandit!). Best single
#                       program: 55.95 on the evolution set, 55.23 on the
#                       held-out 50 (beats the published value on 47 of 50).
#   schedule = :race    successive halving to one survivor, then widening to
#                       lattice neighbours only once the pool has exhausted
#                       its own roots (race_then_widen!). 55.50; wins on about
#                       a fifth of n where :bandit does not.
#
# solve_best(n) runs both and keeps the lower defect -- that is how the
# census values in records/ were produced.
#
# The CSP-call budget is a COUNT, not a wall-clock allowance, so results are
# deterministic: the same n, library and schedule give the same dissection on
# any machine as long as `timeout` does not bind first.


# ------------------------------------------------------------------------------
# 1. Shape schedule
# ------------------------------------------------------------------------------

function recommended_configs(n::Int; deep::Bool=false)
    target = 2.895 * (float(n) ^ 0.472)
    if n < 80
        return [(1, k) for k in 8:14]
    elseif n < 180
        return vcat([(1, k) for k in 11:15], [(2, k) for k in 10:13],
                    [(3, 10), (3, 11), (3, 12)])
    end
    # Measured mean %-below-trend per (R, M) from the 684-n survey. Used
    # only to ORDER the racing pool below -- not to pick a winner.
    quality_table = Dict{Tuple{Int,Int}, Float64}(
        (7,12) => 2.18,  (6,14) => 1.65,  (9,11) => 1.57,
        (5,14) => 1.14,  (5,13) => 1.13,  (8,12) => 1.09,
        (6,13) => 0.90,  (6,11) => 0.88,  (5,11) => 0.71,
        (5,12) => 0.43,  (6,12) => 0.32,  (8,14) => 0.27,
        (10,12) => -0.02, (9,12) => -0.11, (7,14) => -0.14,
        (7,13) => -0.20,  (8,13) => -0.22, (3,14) => -0.25,
        (10,11) => -0.31, (10,13) => -0.33, (4,13) => -0.35,
        (4,12) => -0.50,  (8,11) => -0.53, (9,13) => -0.84,
        (4,11) => -0.87,  (9,14) => -1.65, (3,15) => -1.91,
    )
    function quality(r::Int, m::Int)
        haskey(quality_table, (r, m)) && return quality_table[(r, m)]
        r_pen = -0.28 * ((r - 6.5) / 3.2)^2
        m_center = r <= 6 ? 13.0 : (r <= 8 ? 12.0 : 11.0)
        m_pen = -0.16 * (m - m_center)^2
        return r_pen + m_pen
    end
    lo, hi = 0.55 * target, 1.6 * target
    cands = Tuple{Int,Int,Float64,Float64}[]
    for r in 2:14, m in 8:18
        pieces = r * m
        (pieces < lo || pieces > hi) && continue
        push!(cands, (r, m, quality(r, m), pieces / target))
    end
    isempty(cands) && return [(6, 14), (7, 13), (8, 14), (13, 13)]
    buckets = [(0.55, 0.85, 3), (0.85, 1.15, 4), (1.15, 1.6, 4)]
    ordered = Tuple{Int,Int}[]
    used = Set{Tuple{Int,Int}}()
    for (blo, bhi, npick) in buckets
        in_bucket = [c for c in cands if blo <= c[4] < bhi]
        isempty(in_bucket) && continue
        sort!(in_bucket, by = c -> -c[3])
        for c in in_bucket[1:min(npick, length(in_bucket))]
            rm = (c[1], c[2])
            if !(rm in used)
                push!(ordered, rm); push!(used, rm)
            end
        end
    end
    rest = [c for c in cands if !((c[1], c[2]) in used)]
    sort!(rest, by = c -> -c[3])
    for c in rest
        rm = (c[1], c[2])
        rm in used && continue
        push!(ordered, rm); push!(used, rm)
    end
    keep = deep ? 24 : 16
    return ordered[1:min(keep, length(ordered))]
end

# ------------------------------------------------------------------------------
# 2. Resumable arms: one (R, M) shape, its ranked roots and masks
# ------------------------------------------------------------------------------

mutable struct RaceArm
    r::Int
    k::Int
    weights::Union{Nothing, Vector{Int}}
    ranked_roots::Vector{Topology}
    next_root_idx::Int
    trial_queue::Vector{Tuple{Float64,Int,Vector{IntRect},Vector{Vector{LocalBoxCandidate}}}}
    trial_ptr::Int
    box_cache::Dict{Int, Dict{Tuple{Int,Int}, Vector{LocalBoxCandidate}}}
    target_area::Float64
    best_defect::Float64
    best_proxy::Float64
    calls_used::Int
    exhausted::Bool
    max_roots::Int
    outer_trial_cap::Int
    shortlist_cnt::Int
    keep_per_child::Int
    cands_cap::Int
end
# Rank a root's own outer-boxes by area-spread at the nearest rounding --
# cheap, no CSP call. Shared by uniform arm construction.
function _rank_roots(avail_roots::Vector{Topology}, n::Int, max_roots::Int)
    root_score = Vector{Tuple{Float64, Int}}(undef, length(avail_roots))
    for (ri, rt) in enumerate(avail_roots)
        vv = length(rt.position)
        nearest = 0
        for i in 1:vv
            (rt.position[i] * n - floor(rt.position[i] * n)) >= 0.5 && (nearest |= 1 << (i - 1))
        end
        bx = instantiate_square(rt, n, nearest)
        if bx === nothing
            root_score[ri] = (Inf, ri)
        else
            amin = typemax(Int); amax = typemin(Int)
            for b in bx
                ar = (b[3] - b[1]) * (b[4] - b[2])
                ar < amin && (amin = ar)
                ar > amax && (amax = ar)
            end
            root_score[ri] = (Float64(amax - amin), ri)
        end
    end
    sort!(root_score, by = t -> t[1])
    ranked = [avail_roots[t[2]] for t in root_score if isfinite(t[1])]
    isempty(ranked) && (ranked = avail_roots)
    return ranked[1:min(length(ranked), max_roots)]
end
function make_arm(n::Int, r::Int, k::Int,
                   roots_by_r::Dict{Int, Vector{Topology}}, children_by_k::Dict{Int, Vector{Topology}};
                   max_roots::Int, outer_trial_cap::Int,
                   shortlist_cnt::Int, keep_per_child::Int, cands_cap::Int)
    avail_roots = get(roots_by_r, r, Topology[])
    isempty(avail_roots) && return nothing
    (haskey(children_by_k, k) && !isempty(children_by_k[k])) || return nothing
    target_area = (n * n) / Float64(r * k)
    ranked = _rank_roots(avail_roots, n, max_roots)
    return RaceArm(r, k, nothing, ranked, 1,
                   Tuple{Float64,Int,Vector{IntRect},Vector{Vector{LocalBoxCandidate}}}[], 1,
                   Dict{Int, Dict{Tuple{Int,Int}, Vector{LocalBoxCandidate}}}(),
                   target_area, Inf, Inf, 0, false,
                   max_roots, outer_trial_cap, shortlist_cnt, keep_per_child, cands_cap)
end
# Non-uniform sibling: give the boxes of an already-ranked root unequal
# child counts (summing to the fitted target piece count rather than
# r*k), via equalize_weighted, so box areas become proportional to how
# many children each carries. Wrapped so any failure (missing function,
# unexpected return type, missing child bucket for a given count) is a
# clean no-op -- this can only ever ADD a reachable candidate.
function try_make_nonuniform_arm(n::Int, root_template::Topology, r::Int,
                                  children_by_k::Dict{Int, Vector{Topology}}, target_pieces::Float64,
                                  outer_trial_cap::Int, shortlist_cnt::Int, keep_per_child::Int, cands_cap::Int)
    try
        base = clamp(round(Int, target_pieces / r), 8, 18)
        weights = fill(base, r)
        remainder = round(Int, target_pieces) - base * r
        step = remainder >= 0 ? 1 : -1
        idx = 1
        for _ in 1:abs(remainder)
            weights[idx] = clamp(weights[idx] + step, 8, 18)
            idx = idx % r + 1
        end
        length(unique(weights)) == 1 && return nothing
        for w in weights
            (haskey(children_by_k, w) && !isempty(children_by_k[w])) || return nothing
        end
        new_root = equalize_weighted(root_template, Float64.(weights))
        (new_root isa Topology) || return nothing
        total_pieces = sum(weights)
        target_area = (n * n) / Float64(total_pieces)
        return RaceArm(r, 0, weights, Topology[new_root], 1,
                       Tuple{Float64,Int,Vector{IntRect},Vector{Vector{LocalBoxCandidate}}}[], 1,
                       Dict{Int, Dict{Tuple{Int,Int}, Vector{LocalBoxCandidate}}}(),
                       target_area, Inf, Inf, 0, false,
                       1, outer_trial_cap, shortlist_cnt, keep_per_child, cands_cap)
    catch
        return nothing
    end
end
function lattice_neighbors(r::Int, k::Int)
    cands = Tuple{Int,Int}[]
    for (dr, dk) in ((1,0), (-1,0), (0,1), (0,-1))
        nr, nk = r + dr, k + dk
        (nr >= 2 && 8 <= nk <= 18) && push!(cands, (nr, nk))
    end
    return cands
end
# Lazily (re)builds the ranked trial queue for the next untried root,
# advancing next_root_idx exactly once per attempt so no root is ever
# re-ranked or re-scanned across budget top-ups.
function refill_trials!(arm::RaceArm, n::Int, children_by_k::Dict{Int, Vector{Topology}})
    while arm.next_root_idx <= length(arm.ranked_roots)
        root = arm.ranked_roots[arm.next_root_idx]
        arm.next_root_idx += 1
        masks = balanced_outer_masks(root, n; keep=16)
        queue = Tuple{Float64,Int,Vector{IntRect},Vector{Vector{LocalBoxCandidate}}}[]
        for (box_defect, mask, boxes) in masks
            groups = Vector{Vector{LocalBoxCandidate}}()
            all_viable = true
            for (bi, b) in enumerate(boxes)
                bw, bh = b[3] - b[1], b[4] - b[2]
                k_i = arm.weights === nothing ? arm.k : arm.weights[bi]
                avail = get(children_by_k, k_i, Topology[])
                if isempty(avail)
                    all_viable = false
                    break
                end
                cache_for_k = get!(arm.box_cache, k_i) do
                    Dict{Tuple{Int,Int}, Vector{LocalBoxCandidate}}()
                end
                cands = get_cached_box_candidates(avail, bw, bh, arm.target_area, cache_for_k;
                                                  shortlist_count=arm.shortlist_cnt,
                                                  keep_per_child=arm.keep_per_child,
                                                  cands_cap=arm.cands_cap)
                if isempty(cands)
                    all_viable = false
                    break
                end
                push!(groups, cands)
            end
            all_viable || continue
            approx_defect = maximum(g[1].max_area for g in groups) - minimum(g[1].min_area for g in groups)
            push!(queue, (Float64(approx_defect), mask, boxes, groups))
        end
        sort!(queue, by = x -> x[1])
        queue = queue[1:min(length(queue), arm.outer_trial_cap)]
        if !isempty(queue)
            arm.trial_queue = queue
            arm.trial_ptr = 1
            return true
        end
    end
    return false
end
# Spend up to `budget` CSP calls resuming this arm from wherever it left
# off. Returns (calls_used, perfect_found).
function arm_step!(arm::RaceArm, n::Int, children_by_k::Dict{Int, Vector{Topology}},
                    budget::Int, get_cap::Function, record_solution::Function, method_tag::String)
    calls = 0
    perfect = false
    while calls < budget
        if arm.trial_ptr > length(arm.trial_queue)
            if !refill_trials!(arm, n, children_by_k)
                arm.exhausted = true
                break
            end
        end
        approx_d, mask, boxes, groups = arm.trial_queue[arm.trial_ptr]
        arm.trial_ptr += 1
        arm.best_proxy = min(arm.best_proxy, approx_d)
        current_cap = get_cap()
        if approx_d >= current_cap
            continue
        end
        calls += 1
        arm.calls_used += 1
        tiling = solve_interval_csp(groups, boxes; best_defect_cap=current_cap, time_limit_sec=1.0e9)
        if tiling !== nothing
            d = validate_solution(n, tiling)
            if d !== nothing
                arm.best_defect = min(arm.best_defect, Float64(d))
                record_solution(Solution(n, d, tiling, method_tag))
                if d == 0
                    perfect = true
                    break
                end
            end
        end
    end
    return calls, perfect
end
_arm_rank_key(a::RaceArm) = (isfinite(a.best_defect) ? 0 : 1,
                              isfinite(a.best_defect) ? a.best_defect : a.best_proxy)
function _arm_tag(a::RaceArm)
    a.weights === nothing ? "hierarchical (R=$(a.r), M=$(a.k))" : "hierarchical-nu (R=$(a.r))"
end
# Build fresh arms for the untried lattice neighbours (R+-1, M+-1) of each
# seed arm, plus (for uniform seeds) one non-uniform sibling arm whose box
# child-counts are fitted to the target piece count via equalize_weighted.
# `tried` is mutated in place so no (r,k) pair is ever raced twice.
function try_widen_from!(seed_arms::Vector{RaceArm}, n::Int,
                          roots_by_r::Dict{Int, Vector{Topology}}, children_by_k::Dict{Int, Vector{Topology}},
                          tried::Set{Tuple{Int,Int}}, target_pieces::Float64,
                          max_roots::Int, outer_trial_cap::Int, shortlist_cnt::Int,
                          keep_per_child::Int, cands_cap::Int)
    new_arms = RaceArm[]
    for a in seed_arms
        rep_k = a.weights === nothing ? a.k : round(Int, sum(a.weights) / length(a.weights))
        for (nr, nk) in lattice_neighbors(a.r, rep_k)
            (nr, nk) in tried && continue
            push!(tried, (nr, nk))
            na = make_arm(n, nr, nk, roots_by_r, children_by_k;
                          max_roots=max_roots, outer_trial_cap=outer_trial_cap,
                          shortlist_cnt=shortlist_cnt, keep_per_child=keep_per_child, cands_cap=cands_cap)
            na !== nothing && push!(new_arms, na)
        end
        if a.weights === nothing && !isempty(a.ranked_roots)
            nu = try_make_nonuniform_arm(n, a.ranked_roots[1], a.r, children_by_k, target_pieces,
                                         outer_trial_cap, shortlist_cnt, keep_per_child, cands_cap)
            nu !== nothing && push!(new_arms, nu)
        end
    end
    return new_arms
end

# ------------------------------------------------------------------------------
# 3. Scheduler :bandit
# ------------------------------------------------------------------------------

# Successive halving over a pool of lazy, resumable arms, with widening
# ROLLED INTO the loop instead of run as a separate phase afterward.
# Each round every live arm gets `round_budget` more CSP calls (which
# doubles every round, so total spend per round stays roughly constant as
# the pool halves); ties for arms with no validated tiling yet are broken
# by the best optimistic approx-defect seen (free, already computed while
# deciding admissible pruning). Whenever the surviving frontier for a
# round narrows to <=2 arms (or fully exhausts its own root/mask ladder)
# with budget still unspent, fresh lattice-neighbour arms of the current
# best performer are folded into the NEXT round alongside the survivors --
# so an injected neighbour is itself subject to further halving instead of
# receiving one flat top-up. Returns (final_arms_sorted, calls_used,
# perfect_found).
function run_incremental_bandit!(arms::Vector{RaceArm}, n::Int,
                                  roots_by_r::Dict{Int, Vector{Topology}}, children_by_k::Dict{Int, Vector{Topology}},
                                  budget_limit::Int, get_cap::Function, record_solution::Function,
                                  start_time::Float64, has_timeout::Bool, timeout::Float64,
                                  max_roots::Int, outer_trial_cap::Int, shortlist_cnt::Int,
                                  keep_per_child::Int, cands_cap::Int, target_pieces::Float64,
                                  tried::Set{Tuple{Int,Int}}; max_widen_events::Int=4)
    isempty(arms) && return (arms, 0, false)
    calls_used = 0
    perfect = false
    widen_events = 0
    round_budget = 2
    while !isempty(arms) && calls_used < budget_limit &&
          !(has_timeout && (time() - start_time > timeout))
        for arm in arms
            (has_timeout && (time() - start_time > timeout)) && break
            calls_used >= budget_limit && break
            arm.exhausted && continue
            b = min(round_budget, budget_limit - calls_used)
            b <= 0 && continue
            c, p = arm_step!(arm, n, children_by_k, b, get_cap, record_solution, _arm_tag(arm))
            calls_used += c
            p && (perfect = true)
        end
        perfect && break
        (calls_used >= budget_limit || (has_timeout && (time() - start_time > timeout))) && break
        sort!(arms, by = _arm_rank_key)
        active = [a for a in arms if !a.exhausted]
        if isempty(active)
            widen_events >= max_widen_events && break
            seed = isempty(arms) ? RaceArm[] : [arms[1]]
            new_arms = try_widen_from!(seed, n, roots_by_r, children_by_k, tried, target_pieces,
                                        max_roots, outer_trial_cap, shortlist_cnt, keep_per_child, cands_cap)
            isempty(new_arms) && break
            widen_events += 1
            arms = new_arms
            round_budget *= 2
            continue
        end
        keep_n = max(1, ceil(Int, length(active) / 2))
        survivors = active[1:min(keep_n, length(active))]
        new_arms = RaceArm[]
        if widen_events < max_widen_events && length(survivors) <= 2 && calls_used < budget_limit
            new_arms = try_widen_from!(survivors[1:1], n, roots_by_r, children_by_k, tried, target_pieces,
                                        max_roots, outer_trial_cap, shortlist_cnt, keep_per_child, cands_cap)
            !isempty(new_arms) && (widen_events += 1)
        end
        arms = vcat(survivors, new_arms)
        round_budget *= 2
    end
    sort!(arms, by = _arm_rank_key)
    return (arms, calls_used, perfect)
end

# ------------------------------------------------------------------------------
# 4. Scheduler :race
# ------------------------------------------------------------------------------

# Successive halving over a pool of lazy, resumable arms. Each round keeps
# the per-round total roughly constant while the pool halves, so
# per-survivor budget grows automatically. The final survivor gets
# whatever remains in one commitment. Returns
# (surviving_arms_sorted, calls_used, perfect_found, all_remaining_exhausted).
# The last flag is the widen-on-convergence trigger: true means the race
# stopped because every arm ran out of its own ladder, not because the
# deterministic budget ran out.
function race_pool!(arms::Vector{RaceArm}, n::Int, children_by_k::Dict{Int, Vector{Topology}},
                     budget_limit::Int, get_cap::Function, record_solution::Function,
                     start_time::Float64, has_timeout::Bool, timeout::Float64)
    isempty(arms) && return (arms, 0, false, true)
    calls_used = 0
    perfect = false
    num_rounds = max(1, floor(Int, log2(max(length(arms), 1))) + 1)
    round_no = 0
    while length(arms) > 1
        (has_timeout && (time() - start_time > timeout)) && break
        calls_used >= budget_limit && break
        round_no += 1
        n_i = length(arms)
        rounds_left = max(1, num_rounds - round_no + 1)
        per_round_total = max(n_i, (budget_limit - calls_used) ÷ rounds_left)
        budget_each = max(1, per_round_total ÷ n_i)
        for arm in arms
            (has_timeout && (time() - start_time > timeout)) && break
            calls_used >= budget_limit && break
            arm.exhausted && continue
            b = min(budget_each, budget_limit - calls_used)
            b <= 0 && continue
            c, p = arm_step!(arm, n, children_by_k, b, get_cap, record_solution, _arm_tag(arm))
            calls_used += c
            p && (perfect = true)
        end
        perfect && break
        sort!(arms, by = _arm_rank_key)
        keep_n = max(1, ceil(Int, length(arms) / 2))
        arms = arms[1:min(keep_n, length(arms))]
    end
    if !perfect && !isempty(arms) && calls_used < budget_limit &&
       !(has_timeout && (time() - start_time > timeout))
        sort!(arms, by = _arm_rank_key)
        survivor = arms[1]
        if !survivor.exhausted
            b = budget_limit - calls_used
            c, p = arm_step!(survivor, n, children_by_k, b, get_cap, record_solution, _arm_tag(survivor))
            calls_used += c
            perfect = perfect || p
        end
    end
    sort!(arms, by = _arm_rank_key)
    all_exhausted = all(a -> a.exhausted, arms)
    return (arms, calls_used, perfect, all_exhausted)
end
# Race the pool down to one survivor; if it stopped because every arm ran
# out of its own roots/masks (not because the budget ran out), seed fresh
# lattice-neighbour and non-uniform arms from the top two and race again,
# at most `max_widen_rounds` times. Returns (arms, calls_used, perfect).
function race_then_widen!(arms::Vector{RaceArm}, n::Int,
                          roots_by_r::Dict{Int, Vector{Topology}}, children_by_k::Dict{Int, Vector{Topology}},
                          total_budget::Int, get_cap::Function, record_solution::Function,
                          start_time::Float64, has_timeout::Bool, timeout::Float64,
                          max_roots::Int, outer_trial_cap::Int, shortlist_cnt::Int,
                          keep_per_child::Int, cands_cap::Int, target_pieces::Float64,
                          tried::Set{Tuple{Int,Int}}; max_widen_rounds::Int=4)
    arms, calls_total, perfect_found, all_exhausted = race_pool!(arms, n, children_by_k, total_budget,
                                                                  get_cap, record_solution,
                                                                  start_time, has_timeout, timeout)
    widen_round = 0
    while !perfect_found && calls_total < total_budget && all_exhausted &&
          !(has_timeout && (time() - start_time > timeout)) && widen_round < max_widen_rounds
        widen_round += 1
        seeds = arms[1:min(2, length(arms))]
        new_arms = try_widen_from!(seeds, n, roots_by_r, children_by_k, tried, target_pieces,
                                   max_roots, outer_trial_cap, shortlist_cnt, keep_per_child, cands_cap)
        isempty(new_arms) && break
        arms, used, perfect_found, all_exhausted = race_pool!(new_arms, n, children_by_k,
                                                               total_budget - calls_total,
                                                               get_cap, record_solution,
                                                               start_time, has_timeout, timeout)
        calls_total += used
    end
    return (arms, calls_total, perfect_found)
end

# ------------------------------------------------------------------------------
# 5. Small n (< 180) and the last-resort fallback
# ------------------------------------------------------------------------------

# Rank this root's own outer-boxes by area-spread (cheap, no CSP call), then
# spend up to `call_budget` CSP calls exploring its best roots/masks under
# the given breadth caps. Used only by the small-n (< 180) legacy walk,
# which is a short curated list rather than a wide race, so it is left as
# a plain one-shot sweep. Returns (best_defect_seen, calls_used, perfect).
function attempt_hierarchical_config(n::Int, r::Int, k::Int,
                                      roots_by_r::Dict{Int, Vector{Topology}},
                                      children_by_k::Dict{Int, Vector{Topology}},
                                      call_budget::Int,
                                      max_roots::Int, outer_trial_cap::Int,
                                      shortlist_cnt::Int, keep_per_child::Int, cands_cap::Int,
                                      get_cap::Function, record_solution::Function)
    (haskey(roots_by_r, r) && haskey(children_by_k, k)) || return (Inf, 0, false)
    avail_roots = roots_by_r[r]
    avail_children = children_by_k[k]
    (isempty(avail_roots) || isempty(avail_children) || call_budget <= 0) && return (Inf, 0, false)
    target_area = (n * n) / Float64(r * k)
    box_cache = Dict{Tuple{Int,Int}, Vector{LocalBoxCandidate}}()
    ranked_roots = _rank_roots(avail_roots, n, max_roots)
    calls_used = 0
    best_defect_here = Inf
    for root in ranked_roots
        calls_used >= call_budget && break
        balanced_masks = balanced_outer_masks(root, n; keep=16)
        candidate_outer_masks = Tuple{Int, Int, Int, Vector{IntRect}, Vector{Vector{LocalBoxCandidate}}}[]
        for (box_defect, mask, boxes) in balanced_masks
            groups = Vector{Vector{LocalBoxCandidate}}()
            all_viable = true
            for b in boxes
                bw, bh = b[3] - b[1], b[4] - b[2]
                cands = get_cached_box_candidates(avail_children, bw, bh, target_area, box_cache;
                                                  shortlist_count=shortlist_cnt,
                                                  keep_per_child=keep_per_child,
                                                  cands_cap=cands_cap)
                if isempty(cands)
                    all_viable = false
                    break
                end
                push!(groups, cands)
            end
            all_viable || continue
            approx_defect = maximum(g[1].max_area for g in groups) - minimum(g[1].min_area for g in groups)
            push!(candidate_outer_masks, (approx_defect, box_defect, mask, boxes, groups))
        end
        sort!(candidate_outer_masks, by=x -> (x[1], x[2]))
        trials = candidate_outer_masks[1:min(length(candidate_outer_masks), outer_trial_cap)]
        for (approx_d, _box_d, mask, boxes, groups) in trials
            calls_used >= call_budget && break
            best_defect_here = min(best_defect_here, Float64(approx_d))
            current_cap = get_cap()
            approx_d >= current_cap && continue
            calls_used += 1
            tiling = solve_interval_csp(groups, boxes; best_defect_cap=current_cap, time_limit_sec=1.0e9)
            if tiling !== nothing
                d = validate_solution(n, tiling)
                if d !== nothing
                    best_defect_here = min(best_defect_here, Float64(d))
                    record_solution(Solution(n, d, tiling, "hierarchical (R=$r, M=$k)"))
                    d == 0 && return (0.0, calls_used, true)
                end
            end
        end
    end
    return (best_defect_here, calls_used, false)
end
# Last-resort, budget-free safety net: if the race produced nothing at all
# (rare, but a total miss scores zero for this n), try a handful of
# single-level topologies near the fitted piece-count target.
# screen_topology is cheap and does not consume the CSP-call budget, so
# this cannot compete with or dilute the budget accounted for elsewhere.
function fallback_single_level_sweep!(n::Int, children_by_k::Dict{Int, Vector{Topology}},
                                       record_solution::Function,
                                       start_time::Float64, has_timeout::Bool, timeout::Float64)
    target_k = round(Int, 2.895 * (float(n) ^ 0.472))
    ks = unique(clamp.([target_k - 2, target_k, target_k + 2, 8, 10, 12, 14, 16, 18], 8, 18))
    for k in ks
        haskey(children_by_k, k) || continue
        tops = children_by_k[k]
        for top in tops[1:min(20, length(tops))]
            has_timeout && (time() - start_time > timeout) && return
            defect, _, rects = screen_topology(top, n; max_samples=128)
            if defect !== nothing
                d = validate_solution(n, rects)
                if d !== nothing
                    record_solution(Solution(n, d, rects, "single-level ($(top.id))"))
                end
            end
        end
    end
end

# ------------------------------------------------------------------------------
# 6. Entry points
# ------------------------------------------------------------------------------
function solve(n::Int;
               library::Vector{Floorplan} = Floorplan[],
               outer_roots::Vector{Floorplan} = Floorplan[],
               timeout::Float64 = Inf,
               deep::Bool = false,
               on_update::Union{Function, Nothing} = nothing,
               on_config_start::Union{Function, Nothing} = nothing,
               on_step_progress::Union{Function, Nothing} = nothing,
               schedule::Symbol = :bandit)::Union{Solution, Nothing}
    schedule in (:bandit, :race) || throw(ArgumentError("schedule must be :bandit or :race, got :$schedule"))
    start_time = time()
    has_timeout = isfinite(timeout) && timeout > 0
    raw_children = isempty(library) ? (isfile("data/cache/floorplans.cache") ? load_library("data/cache/floorplans.cache") : Floorplan[]) : library
    root_cache   = isfile("data/cache/outer_roots_expanded.cache") ?
                   "data/cache/outer_roots_expanded.cache" :
                   "data/cache/outer_roots.cache"
    raw_roots    = isempty(outer_roots) ? (isfile(root_cache) ? load_library(root_cache) : Floorplan[]) : outer_roots
    children_by_k = Dict{Int, Vector{Topology}}()
    for fp in raw_children
        t = topology(fp)
        push!(get!(children_by_k, length(t.rects), Topology[]), t)
    end
    roots_by_r = Dict{Int, Vector{Topology}}()
    for fp in raw_roots
        t = topology(fp)
        push!(get!(roots_by_r, length(t.rects), Topology[]), t)
    end
    best_solution = nothing
    function record_solution(sol::Solution)
        if best_solution === nothing || sol.defect < best_solution.defect
            best_solution = sol
            on_update !== nothing && on_update(sol)
        end
    end
    get_cap() = best_solution === nothing ? typemax(Int) : best_solution.defect
    if n < 180
        # Sequential walk over a short curated config list. The two
        # schedules differ only in how an R >= 2 config is explored.
        configs = recommended_configs(n; deep=deep)
        csp_budget_per_config = deep ? 1200 : 65
        csp_budget_total      = deep ? 3600 : 100
        csp_total = 0
        config_time_cap = deep ? 300.0 : 60.0
        for (c_idx, (r, k)) in enumerate(configs)
            has_timeout && (time() - start_time > timeout) && break
            csp_calls = 0
            (csp_total >= csp_budget_total) && break
            best_str = best_solution === nothing ? "none" : string(best_solution.defect)
            on_config_start !== nothing && on_config_start(c_idx, length(configs), r, k, r * k, best_str, round(Int, config_time_cap))
            if r == 1
                haskey(children_by_k, k) || continue
                tops = children_by_k[k]
                for (t_idx, top) in enumerate(tops)
                    has_timeout && (time() - start_time > timeout) && break
                    (csp_calls >= csp_budget_per_config || csp_total >= csp_budget_total) && break
                    csp_calls += 1
                    csp_total += 1
                    defect, _, rects = screen_topology(top, n; max_samples=128)
                    if defect !== nothing
                        d = validate_solution(n, rects)
                        if d !== nothing
                            record_solution(Solution(n, d, rects, "single-level ($(top.id))"))
                            d == 0 && return best_solution
                        end
                    end
                    on_step_progress !== nothing && on_step_progress(t_idx, length(tops), 1, 1, 1, 1, time() - start_time, round(Int, config_time_cap), best_solution === nothing ? typemax(Int) : best_solution.defect)
                end
            else
                budget_left = min(csp_budget_per_config, csp_budget_total - csp_total)
                max_roots       = r >= 3 ? 8 : 10
                outer_trial_cap = r >= 3 ? 4 : 5
                if schedule === :race
                    arm = make_arm(n, r, k, roots_by_r, children_by_k;
                                   max_roots = max_roots, outer_trial_cap = outer_trial_cap,
                                   shortlist_cnt = 36, keep_per_child = 6, cands_cap = 48)
                    if arm !== nothing
                        calls, perfect = arm_step!(arm, n, children_by_k, budget_left, get_cap, record_solution, _arm_tag(arm))
                        csp_total += calls
                        perfect && return best_solution
                    end
                else
                    bdef, used, perfect = attempt_hierarchical_config(n, r, k, roots_by_r, children_by_k,
                                                                        budget_left, max_roots, outer_trial_cap,
                                                                        36, 6, 48, get_cap, record_solution)
                    csp_total += used
                    perfect && return best_solution
                end
            end
        end
        best_solution === nothing && fallback_single_level_sweep!(n, children_by_k, record_solution, start_time, has_timeout, timeout)
        return best_solution
    end
    # n >= 180: successive halving over lazy, resumable arms; the schedule
    # decides how the pool is widened (see sections 3 and 4).
    total_budget    = deep ? 3600 : 100
    max_roots       = deep ? 12 : 9
    outer_trial_cap = deep ? 12 : 6
    shortlist_cnt   = deep ? 80 : 40
    keep_per_child  = deep ? 10 : 7
    cands_cap       = deep ? 80 : 50
    target_pieces   = 2.895 * (float(n) ^ 0.472)
    shape_pool = [c for c in recommended_configs(n; deep=deep) if c[1] >= 2]
    isempty(shape_pool) && (shape_pool = [(6, 14), (7, 13), (8, 14), (9, 11)])
    tried = Set{Tuple{Int,Int}}()
    arms = RaceArm[]
    for (r, k) in shape_pool
        (r, k) in tried && continue
        push!(tried, (r, k))
        a = make_arm(n, r, k, roots_by_r, children_by_k;
                     max_roots=max_roots, outer_trial_cap=outer_trial_cap,
                     shortlist_cnt=shortlist_cnt, keep_per_child=keep_per_child, cands_cap=cands_cap)
        a !== nothing && push!(arms, a)
    end
    if !isempty(arms)
        scheduler! = schedule === :race ? race_then_widen! : run_incremental_bandit!
        _, _, perfect = scheduler!(arms, n, roots_by_r, children_by_k, total_budget,
                                   get_cap, record_solution, start_time, has_timeout, timeout,
                                   max_roots, outer_trial_cap, shortlist_cnt, keep_per_child, cands_cap,
                                   target_pieces, tried)
        perfect && return best_solution
    end
    best_solution === nothing && fallback_single_level_sweep!(n, children_by_k, record_solution, start_time, has_timeout, timeout)
    return best_solution
end

"""
    solve_best(n; kwargs...)

Run `solve` with both schedules and return the lower-defect solution
(`:bandit` on a tie). Each schedule gets the full `timeout`, so this costs
about twice as long as a single `solve`. This is how the census values in
`records/` were produced.
"""
function solve_best(n::Int; kwargs...)::Union{Solution, Nothing}
    best = nothing
    for schedule in (:bandit, :race)
        sol = solve(n; kwargs..., schedule = schedule)
        sol === nothing && continue
        (best === nothing || sol.defect < best.defect) && (best = sol)
    end
    return best
end
