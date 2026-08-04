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
using SHA

# --- Scheme config (points / flux / time axes) ---
include("scheme/SchemeConfig.jl")

# --- FR core ---
include("fr/Points.jl")
include("fr/Correction.jl")
include("fr/Operators.jl")

# --- Mesh / state ---
include("mesh/BoundaryConditions.jl")
include("mesh/Mesh1D.jl")
include("mesh/Mesh2D.jl")
include("mesh/MeshMetrics2D.jl")
include("solvestate/SolutionState.jl")
include("solvestate/SolutionState2D.jl")

# --- Equations ---
include("equations/AbstractEquation.jl")
include("equations/LinearAdvection.jl")
include("equations/Burgers.jl")
include("equations/Euler.jl")
include("equations/LinearAdvection2D.jl")
include("equations/Euler2D.jl")

# --- Fluxes ---
include("flux/Rusanov.jl")
include("flux/HLLC.jl")

# --- Capturing hooks + Persson AV baseline (1D + 2D) ---
include("capturing/Interfaces.jl")
include("capturing/PerssonAV.jl")
include("capturing/PerssonAV2D.jl")

register_method!("persson_av", (; kwargs...) -> PerssonAVMethod(; kwargs...))
include("methods/Registry.jl")

# --- Residual + time ---
include("fr/Residual.jl")
include("fr/Residual2D.jl")
include("time/SSP_RK3.jl")
include("time/SSP_RK2.jl")
include("time/Integrate.jl")

# --- Verification ---
include("verification/schema_keys.jl")
include("verification/Metrics.jl")
include("verification/ExactSod.jl")
include("verification/Scoring.jl")
include("verification/Cases.jl")
include("verification/Cases2D.jl")
include("verification/Report.jl")

# --- Invention loop + experiment log + robustness ---
include("invent/Experiment.jl")
include("invent/Candidate.jl")
include("invent/ExperimentLog.jl")
include("invent/LogAnalytics.jl")
include("invent/Invent.jl")
include("invent/Robustness.jl")
include("invent/Snapshot.jl")

# --- I/O ---
include("io/VTKHighOrder.jl")

# --- CLI ---
include("cli/main.jl")

# Exports
export main_cli
export write_report_skeleton, load_report, validate_report_keys
export DEFAULT_SCORING_WEIGHTS, SCORING_FORMULA_VERSION, SCHEMA_VERSION
export REQUIRED_TOP_LEVEL_KEYS, REQUIRED_SUMMARY_KEYS

export SchemeConfig, DEFAULT_SCHEME, scheme_dict, parse_scheme, time_cfl_guidance

export FROperators, build_operators, n_points
export gauss_legendre_nodes_weights, gauss_lobatto_legendre_nodes_weights, differentiation_matrix
export g_DG_endpoints, g_DG_values_and_derivs, legendre_P

export AbstractBC, PeriodicBC, TransmissiveBC, DirichletBC, ReflectingBC, GhostStateBC
export reflect_conserved, exterior_state
export Mesh1D, physical_coords
export Mesh2D, element_index, element_coords, physical_xy, physical_coords_2d, is_solid
export MeshMetrics2D, build_mesh_metrics, build_mesh_metrics_analytic_wavy
export apply_geometry_warp!, make_wavy_mesh2d, is_curved_mesh, wavy_physical, wavy_partials

export SolutionState, SolutionState2D, allocate_state, set_initial_condition!
export discrete_mass, l2_error

export AbstractEquation, LinearAdvection1D, Burgers1D, Euler1D
export LinearAdvection2D, Euler2D
export physical_flux, physical_flux_x, physical_flux_y, numerical_flux, numerical_flux_n
export max_wave_speed, max_wave_speed_n, n_equations
export rusanov_flux, hllc_flux, hllc_flux_n, interface_flux, interface_flux_n
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
export viscous_mass_residual_scale, viscous_mass_residual_scale_2d
export sensor_field_2d, element_viscosities_2d, element_max_wavespeed_2d

export residual!, ssp_rk3!, ssp_rk3_step!, ssp_rk2!, ssp_rk2_step!, integrate!, compute_dt
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

export run_advection2d_smooth_order, run_euler2d_smooth_order, run_euler2d_discontinuous
export run_m8_2d_suite, run_p31_2d_capturing_suite
export run_freestream_preservation_2d, run_advection2d_curved_order
export run_euler2d_curved_discontinuous, run_p32_curved_suite
export riemann2d_ic, riemann2d_cfg3_ic, RIEMANN2D_CFG3, RIEMANN2D_CFG6, run_euler2d_riemann
export isentropic_vortex_primitives, run_isentropic_vortex_order, run_p33a_benchmark_suite
export double_mach_states, double_mach_ic, run_double_mach_reflection
export make_ffs_solid_mask, run_forward_facing_step, run_p33b_optional_suite

export list_methods, describe_methods, ScaledPerssonMethod
export invent_method, score_reports, run_method_report, write_report
export classify_candidate, print_candidate_summary, score_suite_relative
export tradeoff_ok, DEFAULT_SCORE_MARGIN
export FROZEN_INVENT_SCHEME, default_experiment_log_path, default_experiment_log_yaml_path
export make_entry_id, entry_from_invent, format_entry_markdown
export append_experiment_entry!, invent_append_log!, list_experiment_entry_ids
export narrative_required, package_root, NARRATIVE_PLACEHOLDER
export parse_experiment_log, parse_experiment_log_text, get_experiment_entry
export log_summary, log_frontier, log_pareto, log_lessons
export format_log_summary_text, format_log_frontier_text, format_log_pareto_text
export format_log_lessons_text, format_log_entry_text
export NEAR_MISS_COMPOSITE_MARGIN
export METHOD_SOURCE_MAP, resolve_method_sources, freeze_snapshot, verify_snapshot, snapshot_tables
export SNAPSHOT_SCHEMA_VERSION
export robustness_cells, scheme_slug, evaluate_robustness_cell, run_robustness_matrix
export assess_publication_grade, cell_ok, order_preserved, append_robustness_log_entry!

export write_vtu_high_order, write_vtu_high_order_with_capturing
export compute_capturing_diagnostics_2d, dissipation_operator
export vtk_lagrange_line_nodes, vtk_lagrange_quad_nodes
export vtk_point_counts_1d, vtk_point_counts_2d
export gl_to_equi_interp, parse_vtu_basic, VTK_LAGRANGE_LINE, VTK_LAGRANGE_QUAD

end # module
