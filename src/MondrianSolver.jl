module MondrianSolver

using LinearAlgebra
using Random
using JSON3
using Serialization

# Include all 6 modular components in clean dataflow order
include("core.jl")        # 1. Fundamental data structures & validators
include("io.jl")          # 2. Disk I/O, JSON, OEIS solutions, cache storage
include("geometry.jl")    # 3. Continuous Newton equalizer & atomic mutations
include("generators.jl")  # 4. Planar electrical duals, slicing, & cache expansion
include("csp.jl")         # 5. Discrete rounding & Forward-Checking CSP
include("solver.jl")      # 6. Adaptive strategy & master solve(n) loop

# Core Exports
export Floorplan, Topology, Solution, Rect, IntRect, RectRef, Dim
export LOW, HIGH, canonical_dim
export is_valid, validate_solution, topology, floorplan

# I/O Exports
export floorplan_from_solution, import_floorplans, parse_oeis_quilts
export save_library, load_library, save_solution, update_solutions_json

# Geometry & Mutation Exports
export equalize, equalize_weighted
export t_flips, wall_slides, rectangular_blocks, substitute_block, mutate

# Generator & Expansion Exports
export generate_electrical, generate_electrical_library
export generate_slicing, generate_slicing_library
export expand_electrical_portfolio, generate_quilt_portfolio

# CSP & Solver Exports
export instantiate_square, screen_topology, balanced_outer_masks, solve_interval_csp, solve

end # module MondrianSolver