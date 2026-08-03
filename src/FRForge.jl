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
include("equations/Euler.jl")

# --- Fluxes ---
include("flux/Rusanov.jl")

# --- Capturing hooks + Persson AV baseline ---
include("capturing/Interfaces.jl")
include("capturing/PerssonAV.jl")

# Register Persson after type is defined
register_method!("persson_av", (; kwargs...) -> PerssonAVMethod(; kwargs...))

# Agent methods live only under src/methods/
include("methods/Registry.jl")

# --- Residual + time ---
include("fr/Residual.jl")
include("time/SSP_RK3.jl")

# --- Verification ---
include("verification/schema_keys.jl")
include("verification/Metrics.jl")
include("verification/ExactSod.jl")
include("verification/Scoring.jl")
include("verification/Cases.jl")
include("verification/Report.jl")

# --- Invention loop ---
include("invent/Experiment.jl")
include("invent/Candidate.jl")
include("invent/Invent.jl")

# --- I/O ---
include("io/VTKHighOrder.jl")

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

export AbstractEquation, LinearAdvection1D, Burgers1D, Euler1D
export physical_flux, numerical_flux, max_wave_speed, n_equations
export rusanov_flux
export pressure, velocity, sound_speed, primitives_to_conserved, conserved_to_primitives
export positivity_ok, positivity_ok_state

export AbstractCapturingMethod, NullCapturing
export AbstractShockSensor, AbstractDissipationOperator
export NullSensor, NullDissipation
export PerssonSensor, ElementArtificialViscosity, PerssonAVMethod
export default_persson_params, method_params
export get_capturing_method, register_method!, METHOD_REGISTRY
export sense!, apply_dissipation!, preprocess_state!, extrapolate_interface!
export numerical_flux_method, post_step!
export viscous_mass_residual_scale

export residual!, ssp_rk3!, ssp_rk3_step!, compute_dt
export l2_error_all

export run_advection_smooth_order, run_advection_conservation, run_m1_advection_suite
export run_burgers_conservation, run_burgers_oscillation, run_m2_burgers_suite
export burgers_square_ic
export run_euler_smooth_order, run_euler_conservation, run_m3_euler_suite
export run_bc_transmissive_test, run_bc_dirichlet_test
export euler_density_wave_conserved
export run_persson_vs_null_burgers, run_m4_capturing_suite
export run_sod, run_shu_osher, run_m5_quant_suite, shu_osher_ic, sod_ic
export SodProblem, sod_exact, sod_exact_conserved
export score_suite_absolute, apply_scores!, collect_hard_gate_failures
export sample_solution_1d, shock_thickness_sp, excess_dissipation, smooth_region_mask
export observed_orders, order_pass, solution_extrema, overshoot_metric

export list_methods, describe_methods, ScaledPerssonMethod
export invent_method, score_reports, run_method_report, write_report
export classify_candidate, print_candidate_summary, score_suite_relative
export tradeoff_ok, DEFAULT_SCORE_MARGIN

export write_vtu_high_order, vtk_lagrange_line_nodes, vtk_point_counts_1d
export gl_to_equi_interp, parse_vtu_basic, VTK_LAGRANGE_LINE, VTK_LAGRANGE_QUAD

end # module
