#!/usr/bin/env julia

# ==============================================================================
# plot.jl — High-Resolution Mondrian Solution Visualizer (Luxor.jl + PicoSAT.jl)
# ==============================================================================

using JSON3
using Luxor
using PicoSAT

# Top-level include so methods are defined in the global world age
include(joinpath(@__DIR__, "verify.jl"))

# Strict 4-color Mondrian palette: (fill_hex, text_hex)
const MONDRIAN_PALETTE = [
    ("#F2F1EA", "#111111"),   # 1: Off-white
    ("#D0161B", "#F2F1EA"),   # 2: Red
    ("#1B4C9C", "#F2F1EA"),   # 3: Blue
    ("#F4C700", "#111111"),   # 4: Yellow
]

"""
Build adjacency list: two rectangles are adjacent only if they share a boundary
segment of strictly positive length (touching at a single point does not count).
"""
function build_adjacency(rects::Vector{NTuple{4,Int}})::Vector{Set{Int}}
    n = length(rects)
    adj = [Set{Int}() for _ in 1:n]

    for i in 1:n
        x0i, y0i, x1i, y1i = rects[i]
        for j in i+1:n
            x0j, y0j, x1j, y1j = rects[j]

            # Vertical shared boundary
            if x1i == x0j || x1j == x0i
                if min(y1i, y1j) > max(y0i, y0j)
                    push!(adj[i], j); push!(adj[j], i)
                    continue
                end
            end

            # Horizontal shared boundary
            if y1i == y0j || y1j == y0i
                if min(x1i, x1j) > max(x0i, x0j)
                    push!(adj[i], j); push!(adj[j], i)
                end
            end
        end
    end
    return adj
end

"""
Find an exact 4-coloring via SAT using PicoSAT.jl.
The four-color theorem guarantees a 4-coloring exists for planar contact graphs.
"""
function color_rectangles(adj::Vector{Set{Int}})::Vector{Int}
    n = length(adj)
    vvar(v::Int, c::Int) = (v - 1) * 4 + c

    clauses = Vector{Int}[]

    # 1. Every vertex gets at least one color
    for v in 1:n
        push!(clauses, [vvar(v, 1), vvar(v, 2), vvar(v, 3), vvar(v, 4)])
    end

    # 2. Every vertex gets at most one color
    for v in 1:n, c1 in 1:4, c2 in c1+1:4
        push!(clauses, [-vvar(v, c1), -vvar(v, c2)])
    end

    # 3. Adjacent vertices cannot share the same color
    for u in 1:n
        for v in adj[u]
            if u < v
                for c in 1:4
                    push!(clauses, [-vvar(u, c), -vvar(v, c)])
                end
            end
        end
    end

    sol = PicoSAT.solve(clauses)
    sol isa Symbol && error("PicoSAT could not find a 4-coloring (result: $sol)")

    colors = zeros(Int, n)
    for v in 1:n
        for c in 1:4
            if sol[vvar(v, c)] > 0
                colors[v] = c
                break
            end
        end
    end
    return colors
end

"""Render a high-resolution solution to PNG using Luxor.jl."""
function plot_solution(n::Int, rects::Vector{NTuple{4,Int}}, defect::Int, outpath::String; canvas_size::Int=2400)
    adj = build_adjacency(rects)
    colors = color_rectangles(adj)

    margin_top    = canvas_size * 0.065
    margin_side   = canvas_size * 0.035
    margin_bottom = canvas_size * 0.035
    usable_size   = canvas_size - 2 * margin_side
    scale         = usable_size / n
    canvas_height = round(Int, margin_top + n * scale + margin_bottom)

    Drawing(canvas_size, canvas_height, outpath)
    background("white")

    # Title header
    setcolor("#111111")
    fontsize(canvas_size * 0.022)
    fontface("sans-serif")
    text("n = $n    defect = $defect    pieces = $(length(rects))",
         Point(canvas_size / 2, margin_top / 2), halign=:center, valign=:middle)

    # Resolution-adjusted border stroke
    border_width = max(3.0, min(8.0, 500.0 / n))

    # Render rectangles
    for (i, r) in enumerate(rects)
        x0, y0, x1, y1 = r
        w, h = x1 - x0, y1 - y0
        px = margin_side + x0 * scale
        py = margin_top + y0 * scale
        pw = w * scale
        ph = h * scale

        fill_col, text_col = MONDRIAN_PALETTE[colors[i]]

        # Fill rectangle
        setcolor(fill_col)
        rect(px, py, pw, ph, :fill)

        # Draw black border stroke
        setcolor("#111111")
        setline(border_width)
        rect(px, py, pw, ph, :stroke)

        # High-res label sizing & orientation
        label = "$(w)×$(h)"
        avail_len = max(pw, ph)
        avail_wid = min(pw, ph)

        max_f_len = (avail_len * 0.85) / (0.65 * length(label))
        max_f_wid = avail_wid * 0.70
        max_font_cap = canvas_size * 0.024
        fsize = clamp(min(max_f_len, max_f_wid), 8.0, max_font_cap)

        if fsize >= 10.0
            @layer begin
                translate(Point(px + pw / 2, py + ph / 2))
                h > w && rotate(pi / 2)  # Rotate along long axis
                fontsize(fsize)
                fontface("sans-serif")
                setcolor(text_col)
                text(label, O, halign=:center, valign=:middle)
            end
        end
    end

    finish()
end

function parse_args(args::Vector{String})
    file = "solutions.json"
    outdir = "plots"
    canvas_size = 2400
    targets = String[]

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("-h", "--help")
            println("Usage: julia plot.jl [FILE] [--outdir DIR] [--n N1,N2,...] [--size PIXELS]")
            exit(0)
        elseif arg == "--outdir"
            i += 1; outdir = args[i]
        elseif arg == "--size"
            i += 1; canvas_size = parse(Int, args[i])
        elseif arg == "--n"
            i += 1; append!(targets, split(args[i], ","))
        elseif !startswith(arg, "-")
            file = arg
        end
        i += 1
    end
    return file, outdir, targets, canvas_size
end

function main()
    file, outdir, targets, canvas_size = parse_args(ARGS)

    if !isfile(file)
        println(stderr, "Error: $file not found.")
        return 1
    end

    data = JSON3.read(read(file, String))
    mkpath(outdir)

    keys_to_plot = isempty(targets) || "all" in targets ?
        sort(collect(keys(data)), by=k -> parse(Int, String(k))) :
        filter(k -> String(k) in targets, collect(keys(data)))

    for k in keys_to_plot
        n = parse(Int, String(k))
        defect, res = verify_solution(n, data[k])

        if defect === nothing
            println("Error a($n): $res (skipping plot)")
            continue
        end

        outpath = joinpath(outdir, "mondrian_$n.png")
        plot_solution(n, res, defect, outpath; canvas_size=canvas_size)
        println("a($n)=$defect  ->  saved $outpath ($(canvas_size)px)")
    end
    return 0
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end