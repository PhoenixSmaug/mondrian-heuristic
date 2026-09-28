module MondrianSolver

using LinearAlgebra
using Random
using JSON3
using Serialization

# Include all 7 modular components in clean dataflow order
include("core.jl")        # 1. Fundamental data structures & validators
include("io.jl")          # 2. Disk I/O, JSON, OEIS solutions, cache storage
include("geometry.jl")    # 3. Continuous Newton equalizer & atomic mutations
include("prime.jl")       # 4. A340984 prime rectangulations: test & generation
include("generators.jl")  # 5. Planar electrical duals, slicing, & cache expansion
include("csp.jl")         # 6. Discrete rounding & Forward-Checking CSP
include("solver.jl")      # 7. Adaptive strategy & master solve(n) loop

# Core Exports
export Floorplan, Topology, Solution, Rect, IntRect, RectRef, Dim
export LOW, HIGH, canonical_dim
export is_valid, validate_solution, topology, floorplan

# I/O Exports
export floorplan_from_solution, floorplan_from_box, import_floorplans
export parse_oeis_quilts, parse_a340984
export save_library, load_library, save_solution, update_solutions_json
export save_a340984_jsonl

# Geometry & Mutation Exports
export equalize, equalize_weighted
export t_flips, wall_slides, rectangular_blocks, substitute_block, mutate
export dihedral_transform

# Prime (A340984) Exports
export has_rect_block, is_prime_rectangulation, d4_signature
export split_rectangle, grow_prime, generate_prime_level, generate_prime_ladder
export grid_drawing

# Generator & Expansion Exports
export generate_electrical, generate_electrical_library
export generate_slicing, generate_slicing_library
export expand_portfolio_by_k, expand_electrical_portfolio
export generate_quilt_portfolio, generate_prime_portfolio, interleave

# CSP & Solver Exports
export instantiate_square, screen_topology, balanced_outer_masks, solve_interval_csp, solve

end # module MondrianSolver