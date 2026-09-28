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
    floorplan_from_box(raw_rects, w, h; id="unnamed") -> Floorplan

Normalize integer rectangles tiling a `w × h` box into a unit `Floorplan`.
Unlike `floorplan_from_solution` the two axes are scaled independently, so
non-square source boxes are supported.
"""
function floorplan_from_box(raw_rects, w::Real, h::Real; id::String="unnamed")::Floorplan
    sx, sy = Float64(w), Float64(h)
    (sx > 0 && sy > 0) || error("degenerate box: $(w) × $(h)")
    rects = [(Float64(r[1]) / sx, Float64(r[2]) / sy,
              Float64(r[3]) / sx, Float64(r[4]) / sy) for r in raw_rects]
    sort!(rects)
    return Floorplan(id, rects)
end

"""
    reservoir_sample_lines(path, cap, rng; skip) -> Vector{String}

Uniformly sample at most `cap` non-skipped lines of `path` in a single streaming
pass, so that multi-gigabyte enumeration files never have to be held in memory.
"""
function reservoir_sample_lines(path::String, cap::Int, rng::AbstractRNG;
                                skip::Function = _ -> false)::Vector{String}
    cap <= 0 && return String[]
    kept = String[]
    seen = 0
    for line in eachline(path)
        (isempty(line) || skip(line)) && continue
        seen += 1
        if length(kept) < cap
            push!(kept, line)
        else
            j = rand(rng, 1:seen)
            j <= cap && (kept[j] = line)
        end
    end
    return kept
end

"""
    parse_a340984(path; k_range=5:18, max_per_k=typemax(Int), include_crossed=false, rng)

Ingest the prime rectangulations catalogued by OEIS
[A340984](https://oeis.org/A340984) as unit `Floorplan`s.

`path` is either a directory holding `representatives_nK.jsonl` /
`generated_nK.jsonl` files or a single such file. Every record is one D4
equivalence class of a tiling in which no proper sub-collection of 2 .. k-1
rectangles forms a rectangle; in particular no two rectangles share a complete
face. These are exactly the topologies that cannot be assembled by gluing
smaller rectangulations, which makes them a natural complement to the slicing,
electrical and quilt seeds.

`max_per_k` caps (by uniform reservoir sampling) how many classes are taken per
piece count, since a(13) = 16143 and a(14) = 88112; pass `0` for no cap. Records
flagged `crossed` (four rectangles meeting at a point) have fewer than `k - 1`
wall segments and therefore admit no equal-area realization; they are dropped
unless `include_crossed` is set.
"""
function parse_a340984(path::String;
                       k_range::AbstractUnitRange{Int} = 5:18,
                       max_per_k::Int = typemax(Int),
                       include_crossed::Bool = false,
                       rng::AbstractRNG = Random.default_rng())::Vector{Floorplan}
    max_per_k <= 0 && (max_per_k = typemax(Int))
    files = String[]
    if isdir(path)
        for name in sort(readdir(path))
            m = match(r"^(?:representatives|generated)_n(\d+)\.jsonl$", name)
            m === nothing && continue
            parse(Int, m.captures[1]) in k_range && push!(files, joinpath(path, name))
        end
        isempty(files) && return Floorplan[]
    elseif isfile(path)
        push!(files, path)
    else
        error("A340984 source not found: $path")
    end

    # Cheap textual pre-filter for the reservoir; the parsed record is re-checked
    # below so a reformatted source can never smuggle a crossed tiling through.
    skip_crossed = include_crossed ? (_ -> false) : (l -> occursin("\"crossed\":true", l))

    floorplans = Floorplan[]
    for file in files
        for line in reservoir_sample_lines(file, max_per_k, rng; skip=skip_crossed)
            rec = try
                JSON3.read(line)
            catch
                continue
            end
            (!include_crossed && get(rec, :crossed, false)) && continue
            (haskey(rec, :rectangles) && haskey(rec, :width) && haskey(rec, :height)) || continue

            raw = rec.rectangles
            length(raw) in k_range || continue
            id = haskey(rec, :class_id) ? "prime_$(rec.class_id)" :
                 "prime_k$(length(raw))_$(length(floorplans) + 1)"

            fp = floorplan_from_box(raw, rec.width, rec.height; id=id)
            is_valid(fp) && push!(floorplans, fp)
        end
    end

    sort!(floorplans, by=fp -> (length(fp), fp.id))
    return floorplans
end


"""
    save_a340984_jsonl(dir, library; prefix="generated") -> Vector{String}

Write `library` back out in the A340984 representative schema, one
`<prefix>_nK.jsonl` file per piece count, using each floorplan's grid drawing.
The files are readable by the catalogue's own `render.py`. Returns the paths
written.
"""
function save_a340984_jsonl(dir::String, library::Vector{Floorplan};
                            prefix::String="generated")::Vector{String}
    mkpath(dir)
    by_k = Dict{Int, Vector{Floorplan}}()
    for fp in library
        push!(get!(by_k, length(fp), Floorplan[]), fp)
    end

    written = String[]
    for k in sort(collect(keys(by_k)))
        path = joinpath(dir, "$(prefix)_n$(k).jsonl")
        open(path, "w") do io
            for (i, fp) in enumerate(by_k[k])
                w, h, drawn = grid_drawing(fp)
                JSON3.write(io, (schema = "A340984-generated-v1",
                                 class_id = "GEN-$k-$(lpad(i, 6, '0'))",
                                 n = k,
                                 crossed = false,
                                 width = w,
                                 height = h,
                                 rectangles = [[r[1], r[2], r[3], r[4]] for r in drawn],
                                 source = (generator = "mondrian.jl gen-prime",
                                           method = "split-and-repair",
                                           id = fp.id)))
                println(io)
            end
        end
        push!(written, path)
    end
    return written
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