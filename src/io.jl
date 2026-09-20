# ==============================================================================
# io.jl — All File Ingestion, JSON Parsing, OEIS Solutions, & Cache Storage
# ==============================================================================

"""
    floorplan_from_solution(raw_rects, [n]; id="unnamed") -> Floorplan

Normalize integer rectangles from an `n × n` partition into a unit `Floorplan`.
"""
function floorplan_from_solution(raw_rects, n::Union{Real, Nothing}=nothing; id::String="unnamed")
    rects = [ntuple(j -> Float64(r[j]), 4) for r in raw_rects]
    scale = if n !== nothing
        Float64(n)
    else
        max_coord = maximum(max(r[3], r[4]) for r in rects)
        max_coord <= 1.0 ? 1.0 : max_coord
    end
    normalized = [(r[1] / scale, r[2] / scale, r[3] / scale, r[4] / scale) for r in rects]
    sort!(normalized)
    return Floorplan(id, normalized)
end

"""
    import_floorplans(path::String) -> Vector{Floorplan}

Import floorplans from JSON solution files (e.g. `ed_pegg_jr.json` or `solutions.json`).
"""
function import_floorplans(path::String)::Vector{Floorplan}
    isfile(path) || error("file not found: $path")
    data = JSON3.read(read(path, String))
    floorplans = Floorplan[]

    if haskey(data, :rectangles) || haskey(data, "rectangles")
        raw = haskey(data, :rectangles) ? data[:rectangles] : data["rectangles"]
        n = haskey(data, :n) ? data[:n] : (haskey(data, "n") ? data["n"] : nothing)
        id = n !== nothing ? "n=$n" : splitext(basename(path))[1]
        push!(floorplans, floorplan_from_solution(raw, n; id=String(id)))
    else
        for (k, raw) in pairs(data)
            n = tryparse(Int, String(k))
            id = n !== nothing ? "n=$n" : String(k)
            push!(floorplans, floorplan_from_solution(raw, n; id))
        end
    end

    sort!(floorplans, by=fp -> (length(fp), fp.id))
    return floorplans
end

"""
    parse_oeis_quilts(path; max_n=96) -> Vector{Floorplan}

Parse ASCII rectangulation quilts from OEIS `oeis_solutions.txt` (orders 3 through `max_n`).
Extracts contiguous rectangular character blocks and normalizes them to `[0, 1]²`.
"""
function parse_oeis_quilts(path::String; max_n::Int=96)::Vector{Floorplan}
    isfile(path) || error("OEIS solutions file not found: $path")
    lines = readlines(path)
    floorplans = Floorplan[]
    i = 1

    while i <= length(lines)
        line = strip(lines[i])
        m = match(r"^a\((\d+)\)(?:=|<=)(\d+)$", line)
        if m === nothing
            i += 1
            continue
        end

        n = parse(Int, m.captures[1])
        if n > max_n || i + n > length(lines)
            i += 1
            continue
        end

        rows = lines[i+1:i+n]
        labels = sort!(unique(collect(join(rows))))
        rects = Rect[]
        is_solid = true

        for label in labels
            xs, ys = Int[], Int[]
            for y in 1:n, x in 1:n
                if rows[y][x] == label
                    push!(xs, x)
                    push!(ys, y)
                end
            end
            isempty(xs) && continue
            x0, x1 = minimum(xs) - 1, maximum(xs)
            y0, y1 = minimum(ys) - 1, maximum(ys)

            if length(xs) != (x1 - x0) * (y1 - y0)
                is_solid = false
                break
            end
            push!(rects, (x0 / n, y0 / n, x1 / n, y1 / n))
        end

        if is_solid && length(rects) >= 2
            sort!(rects)
            push!(floorplans, Floorplan("oeis_n$(n)_k$(length(rects))", rects))
        end
        i += n + 1
    end

    return floorplans
end

"""
    save_library(path::String, library::Vector{Floorplan})

Save a floorplan library to JSON (`.json`) or binary cache (`.cache` / `.jls`).
"""
function save_library(path::String, library::Vector{Floorplan})::String
    mkpath(dirname(abspath(path)))
    if endswith(lowercase(path), ".json")
        dict = Dict(fp.id => [[r[1], r[2], r[3], r[4]] for r in fp.rectangles] for fp in library)
        open(path, "w") do io
            JSON3.pretty(io, dict)
        end
    else
        serialize(path, library)
    end
    return path
end

"""
    load_library(path::String) -> Vector{Floorplan}

Load a floorplan library from a JSON (`.json`) or binary cache file.
"""
function load_library(path::String)::Vector{Floorplan}
    isfile(path) || error("library file not found: $path")
    if endswith(lowercase(path), ".json")
        return import_floorplans(path)
    else
        data = deserialize(path)
        data isa Vector{Floorplan} || error("file does not contain a Vector{Floorplan}")
        return data
    end
end

"""
    update_solutions_json(path::String, sol::Solution) -> (Bool, Int)

Update the consolidated `solutions.json` file with `sol`.
Only overwrites an existing solution for `sol.n` if `sol.defect` is strictly better (smaller).
Returns `(updated::Bool, current_best_defect::Int)`.
"""
function update_solutions_json(path::String, sol::Solution)::Tuple{Bool, Int}
    mkpath(dirname(abspath(path)))
    data = Dict{String, Any}()

    if isfile(path)
        try
            raw = JSON3.read(read(path, String))
            for (k, v) in pairs(raw)
                data[String(k)] = v
            end
        catch
        end
    end

    key = string(sol.n)
    if haskey(data, key)
        existing = data[key]
        ex_rects = existing isa AbstractVector ? existing : existing["rectangles"]
        ex_areas = [(Int(r[3]) - Int(r[1])) * (Int(r[4]) - Int(r[2])) for r in ex_rects]
        ex_defect = maximum(ex_areas) - minimum(ex_areas)
        if sol.defect >= ex_defect
            return (false, ex_defect)
        end
    end

    data[key] = [[r[1], r[2], r[3], r[4]] for r in sol.rectangles]
    open(path, "w") do io
        JSON3.pretty(io, data)
        println(io)
    end
    return (true, sol.defect)
end

"""
    save_solution(path::String, sol::Solution)

Export an integer solution directly to `path` or update `solutions.json`.
"""
function save_solution(path::String, sol::Solution)::String
    if basename(path) == "solutions.json"
        update_solutions_json(path, sol)
        return path
    end
    mkpath(dirname(abspath(path)))
    payload = Dict{String, Any}(
        "n" => sol.n,
        "defect" => sol.defect,
        "method" => sol.method,
        "rectangles" => [[r[1], r[2], r[3], r[4]] for r in sol.rectangles]
    )
    open(path, "w") do io
        JSON3.pretty(io, payload)
    end
    return path
end