"""
    FRForge

High-order Flux Reconstruction laboratory for inventing and evaluating
novel shock-capturing methods on the compressible Euler equations.

This package is intentionally green-field: FR operators, residual evaluation,
verification, and I/O are implemented from first principles with clarity and
verifiability prioritized over performance in early milestones.
"""
module FRForge

using LinearAlgebra
using ArgParse
using Dates
using JSON

# --- FR core ---
include("fr/Points.jl")
include("fr/Correction.jl")
include("fr/Operators.jl")

# --- Mesh / state ---
include("mesh/BoundaryConditions.jl")
include("mesh/Mesh1D.jl")
include("solvestate/SolutionState.jl")

# --- Equations ---
include("equations/AbstractEquation.jl")
include("equations/LinearAdvection.jl")
include("equations/Burgers.jl")

# --- Fluxes ---
include("flux/Rusanov.jl")

# --- Capturing hooks (nulls) ---
include("capturing/Interfaces.jl")

# --- Residual + time ---
include("fr/Residual.jl")
include("time/SSP_RK3.jl")

# --- Verification ---
include("verification/schema_keys.jl")
include("verification/Metrics.jl")
include("verification/Cases.jl")
include("verification/Report.jl")

# --- CLI ---
include("cli/main.jl")

# Exports
export main_cli
export write_report_skeleton, load_report, validate_report_keys
export DEFAULT_SCORING_WEIGHTS, SCORING_FORMULA_VERSION, SCHEMA_VERSION
export REQUIRED_TOP_LEVEL_KEYS, REQUIRED_SUMMARY_KEYS

export FROperators, build_operators, n_points
export gauss_legendre_nodes_weights, differentiation_matrix
export g_DG_endpoints, g_DG_values_and_derivs, legendre_P

export AbstractBC, PeriodicBC, TransmissiveBC, DirichletBC
export Mesh1D, physical_coords

export SolutionState, allocate_state, set_initial_condition!
export discrete_mass, l2_error

export AbstractEquation, LinearAdvection1D, Burgers1D
export physical_flux, numerical_flux, max_wave_speed, n_equations
export rusanov_flux

export AbstractCapturingMethod, NullCapturing
export AbstractShockSensor, AbstractDissipationOperator
export NullSensor, NullDissipation

export residual!, ssp_rk3!, ssp_rk3_step!, compute_dt

export run_advection_smooth_order, run_advection_conservation, run_m1_advection_suite
export run_burgers_conservation, run_burgers_oscillation, run_m2_burgers_suite
export burgers_square_ic
export observed_orders, order_pass, solution_extrema, overshoot_metric

end # module
