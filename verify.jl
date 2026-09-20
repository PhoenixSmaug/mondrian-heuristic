#!/usr/bin/env julia

# ==============================================================================
# verify.jl — Fast Memory-Independent Mondrian Solution Validator
# ==============================================================================

using JSON3

function verify_solution(n::Int, raw_rects)
    rects = raw_rects isa AbstractVector ? raw_rects : raw_rects["rectangles"]
    length(rects) >= 2 || return (nothing, "Trivial solution (only one rectangle)")

    seen_dims = Set{Tuple{Int,Int}}()
    total_area = 0
    min_a, max_a = typemax(Int), 0
    parsed = NTuple{4,Int}[]

    # 1. Check coordinates, positive dimensions, and pairwise non-congruence
    for (idx, r) in enumerate(rects)
        x0, y0, x1, y1 = Int(r[1]), Int(r[2]), Int(r[3]), Int(r[4])
        (0 <= x0 < x1 <= n && 0 <= y0 < y1 <= n) ||
            return (nothing, "Rectangle $idx out of bounds: [$x0, $y0, $x1, $y1]")

        w, h = x1 - x0, y1 - y0
        dim_pair = minmax(w, h)
        dim_pair in seen_dims && return (nothing, "Duplicate rectangle dimensions $(w)x$(h)")
        push!(seen_dims, dim_pair)

        a = w * h
        total_area += a
        min_a = min(min_a, a)
        max_a = max(max_a, a)
        push!(parsed, (x0, y0, x1, y1))
    end

    # 2. Pairwise bounding-box overlap check (no n x n grid allocation)
    for i in 1:length(parsed), j in 1:i-1
        r1, r2 = parsed[i], parsed[j]
        if max(r1[1], r2[1]) < min(r1[3], r2[3]) && max(r1[2], r2[2]) < min(r1[4], r2[4])
            return (nothing, "Overlapping rectangles at index $j and $i")
        end
    end

    # 3. Exact area coverage
    total_area == n * n || return (nothing, "Incomplete tiling (area sum $total_area != $(n*n))")

    return (max_a - min_a, parsed)
end

function main()
    file = get(ARGS, 1, "solutions.json")
    if !isfile(file)
        println(stderr, "Error: $file not found.")
        return 1
    end

    data = JSON3.read(read(file, String))
    sorted_keys = sort(collect(keys(data)), by=k -> parse(Int, String(k)))

    for k in sorted_keys
        n = parse(Int, String(k))
        defect, err = verify_solution(n, data[k])
        if defect !== nothing
            println("M($n)=$defect")
        else
            println("Error M($n): $err")
        end
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end