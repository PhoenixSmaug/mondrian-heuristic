# ==============================================================================
# solver.jl — Adaptive Strategy & Master Solve Coordination
# ==============================================================================

function recommended_configs(n::Int; deep::Bool=false)
    if n < 80
        return [(1, k) for k in 8:14]
    elseif n < 180
        return vcat([(1, k) for k in 11:15], [(2, k) for k in 10:13], [(3, 10), (3, 11), (3, 12)])
    elseif n < 350
        base = [(3, 11), (3, 12), (3, 13), (4, 11), (4, 12), (4, 13), (5, 10), (5, 11)]
        return deep ? vcat(base, [(3, 14), (4, 14), (5, 12)]) : base
    elseif n < 700
        # High piece counts prioritized: R=6,M=13 (78), R=7,M=11 (77), R=13,M=10 (130), R=13,M=13 (169)
        base = [(6, 13), (5, 14), (7, 11), (8, 14), (13, 10), (13, 13)]
        return deep ? vcat(base, [(7, 13), (6, 14), (8, 13), (14, 14), (5, 13), (7, 12)]) : base
    else
        base = [(6, 14), (7, 13), (8, 14), (13, 13), (14, 14)]
        return deep ? vcat(base, [(13, 14), (14, 16), (9, 14)]) : base
    end
end

"""
    solve(n; [library], [outer_roots], [timeout=Inf], [deep=false], [on_update], [on_config_start], [on_step_progress])

Master solver with a generous 5-minute budget per configuration, live progress bars, and immediate saving.
"""
function solve(n::Int;
               library::Vector{Floorplan} = Floorplan[],
               outer_roots::Vector{Floorplan} = Floorplan[],
               timeout::Float64 = Inf,
               deep::Bool = false,
               on_update::Union{Function, Nothing} = nothing,
               on_config_start::Union{Function, Nothing} = nothing,
               on_step_progress::Union{Function, Nothing} = nothing)::Union{Solution, Nothing}
    start_time = time()
    has_timeout = isfinite(timeout) && timeout > 0

    raw_children = isempty(library) ? (isfile("data/cache/floorplans.cache") ? load_library("data/cache/floorplans.cache") : Floorplan[]) : library
    raw_roots    = isempty(outer_roots) ? (isfile("data/cache/outer_roots.cache") ? load_library("data/cache/outer_roots.cache") : Floorplan[]) : outer_roots

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

    configs = recommended_configs(n; deep=deep)

    # 5 minutes per config in deep mode (300s), 1 minute in standard mode (60s)
    config_time_cap = deep ? 300.0 : 60.0

    for (c_idx, (r, k)) in enumerate(configs)
        has_timeout && (time() - start_time > timeout) && break
        config_start = time()

        best_str = best_solution === nothing ? "none" : string(best_solution.defect)
        on_config_start !== nothing && on_config_start(c_idx, length(configs), r, k, r * k, best_str, round(Int, config_time_cap))

        # ── 1. Single-Level (R = 1) ──
        if r == 1
            haskey(children_by_k, k) || continue
            tops = children_by_k[k]
            for (t_idx, top) in enumerate(tops)
                has_timeout && (time() - start_time > timeout) && break
                (time() - config_start > config_time_cap) && break

                defect, _, rects = screen_topology(top, n; max_samples=128)
                if defect !== nothing
                    d = validate_solution(n, rects)
                    if d !== nothing
                        record_solution(Solution(n, d, rects, "single-level ($(top.id))"))
                        d == 0 && return best_solution
                    end
                end
                on_step_progress !== nothing && on_step_progress(t_idx, length(tops), 1, 1, 1, 1, time() - config_start, round(Int, config_time_cap), best_solution === nothing ? typemax(Int) : best_solution.defect)
            end
        else
            # ── 2. Hierarchical (R >= 2) ──
            haskey(roots_by_r, r) || continue
            haskey(children_by_k, k) || continue
            avail_roots = roots_by_r[r]
            avail_children = children_by_k[k]
            target_area = (n * n) / Float64(r * k)

            box_cache = Dict{Tuple{Int,Int}, Vector{LocalBoxCandidate}}()

            # Scaled roots and trial depth for 5-minute budget
            max_roots       = r >= 11 ? (deep ? 6 : 3)  : (r >= 7 ? (deep ? 10 : 5) : (deep ? 12 : 6))
            outer_trial_cap = r >= 11 ? (deep ? 6 : 3)  : (r >= 7 ? (deep ? 12 : 6) : (deep ? 16 : 8))
            shortlist_cnt   = deep ? 80 : 36
            keep_per_child  = deep ? 10 : 6
            cands_cap       = deep ? 80 : 48
            time_limit_csp  = deep ? 10.0 : 4.0

            roots_to_test = avail_roots[1:min(length(avail_roots), max_roots)]
            total_steps = length(roots_to_test) * outer_trial_cap
            current_step = 0

            for (root_idx, root) in enumerate(roots_to_test)
                has_timeout && (time() - start_time > timeout) && break
                (time() - config_start > config_time_cap) && break

                balanced_masks = balanced_outer_masks(root, n; keep=(deep ? 24 : 16))
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

                for (trial_idx, (approx_d, _, mask, boxes, groups)) in enumerate(trials)
                    has_timeout && (time() - start_time > timeout) && break
                    (time() - config_start > config_time_cap) && break
                    current_step += 1

                    current_cap = best_solution === nothing ? typemax(Int) : best_solution.defect

                    on_step_progress !== nothing && on_step_progress(
                        current_step, total_steps,
                        root_idx, length(roots_to_test),
                        trial_idx, length(trials),
                        time() - config_start,
                        round(Int, config_time_cap),
                        current_cap
                    )

                    approx_d >= current_cap && continue

                    tiling = solve_interval_csp(groups, boxes; best_defect_cap=current_cap, time_limit_sec=time_limit_csp)
                    if tiling !== nothing
                        d = validate_solution(n, tiling)
                        if d !== nothing
                            record_solution(Solution(n, d, tiling, "hierarchical (R=$r, M=$k)"))
                            d == 0 && return best_solution
                        end
                    end
                end
            end
        end
    end

    return best_solution
end