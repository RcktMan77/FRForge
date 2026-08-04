"""
    FRForge

High-order Flux Reconstruction laboratory for inventing and evaluating
novel shock-capturing methods on the compressible Euler equations.

FR operators, residual evaluation, verification, and I/O are implemented from
first principles. The default invent residual path is serial and bit-deterministic;
optional multi-thread residual kernels are docs/long-run only.
"""
module FRForge

using LinearAlgebra
using Base.Threads
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
include("fr/Threading.jl")

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
include("fr/ResidualWorkspace.jl")
include("fr/Residual.jl")
include("fr/Residual2D.jl")
include("time/SSP_RK3.jl")
include("time/SSP_RK2.jl")
include("time/Integrate.jl")

# --- Verification ---
include("verification/schema_keys.jl")
include("verification/Metrics.jl")
include("verification/CaseReport.jl")
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

# Confirm after VTK (optional --vtk path uses HO writer)
include("invent/Confirm.jl")

# --- CLI (constants → commands → dispatcher) ---
include("cli/scheme_args.jl")
include("cli/constants.jl")
include("cli/run.jl")
include("cli/test_cmd.jl")
include("cli/invent_cmd.jl")
include("cli/log_cmd.jl")
include("cli/snapshot_cmd.jl")
include("cli/main.jl")

# =============================================================================
# Export surface
#
# Convention:
#   - Names without a leading underscore are the supported `using FRForge` API
#     (CLI, invent agents, verification, and unit tests).
#   - Leading underscore (`_foo`) means package-internal; use `FRForge._foo` only
#     from tests when necessary. Prefer not to depend on them.
#   - Sections below group the public surface by audience.
# =============================================================================

# --- Process / reports ---
export main_cli
export write_report_skeleton, load_report, validate_report_keys, write_report
export write_json_pretty, report_trio_paths, stamp_workflow_report!, print_report_trio
export resolve_log_paths, report_artifact_dict, maybe_append_workflow_log!
export DEFAULT_SCORING_WEIGHTS, SCORING_FORMULA_VERSION, SCHEMA_VERSION
export REQUIRED_TOP_LEVEL_KEYS, REQUIRED_SUMMARY_KEYS
export case_report_dict

# --- Scheme axes (invent freezes DEFAULT_SCHEME) ---
export SchemeConfig, DEFAULT_SCHEME, scheme_dict, parse_scheme, time_cfl_guidance

# --- FR operators & residual time marching ---
export FROperators, build_operators, n_points
export gauss_legendre_nodes_weights, gauss_lobatto_legendre_nodes_weights, differentiation_matrix
export g_DG_endpoints, g_DG_values_and_derivs, legendre_P
export residual!, ssp_rk3!, ssp_rk3_step!, ssp_rk2!, ssp_rk2_step!, integrate!, compute_dt
export frforge_thread_count, residual_threading_enabled, with_frforge_threads, with_serial_residual
# Threading loop helper (tests / advanced residual work)
export foreach_element

# --- Mesh / BC / state ---
export AbstractBC, PeriodicBC, TransmissiveBC, DirichletBC, ReflectingBC, GhostStateBC
export reflect_conserved, exterior_state
export Mesh1D, physical_coords
export Mesh2D, element_index, element_coords, physical_xy, physical_coords_2d, is_solid
export MeshMetrics2D, build_mesh_metrics, build_mesh_metrics_analytic_wavy
export apply_geometry_warp!, make_wavy_mesh2d, is_curved_mesh
export SolutionState, SolutionState2D, allocate_state, set_initial_condition!
export discrete_mass, l2_error
# l2_error_all is internal (FRForge.l2_error_all); not exported

# --- Equations & fluxes ---
export AbstractEquation, LinearAdvection1D, Burgers1D, Euler1D
export LinearAdvection2D, Euler2D
export physical_flux, physical_flux!, physical_flux_x, physical_flux_y, numerical_flux, numerical_flux_n
export max_wave_speed, max_wave_speed_n, n_equations
export rusanov_flux, hllc_flux, hllc_flux_n, interface_flux, interface_flux_n
export pressure, velocity, sound_speed, primitives_to_conserved, conserved_to_primitives
export positivity_ok, positivity_ok_state

# --- Capturing hooks & registry (primary invent surface) ---
export AbstractCapturingMethod, NullCapturing
export AbstractShockSensor, AbstractDissipationOperator
export NullSensor, NullDissipation
export PerssonSensor, ElementArtificialViscosity, PerssonAVMethod
export default_persson_params, method_params
export get_capturing_method,
    register_method!, require_registered_method, list_methods, describe_methods
export sense!, apply_dissipation!, preprocess_state!, extrapolate_interface!
export numerical_flux_method, post_step!
export ScaledPerssonMethod
export viscous_mass_residual_scale, viscous_mass_residual_scale_2d
# METHOD_REGISTRY is internal; use list_methods / get_capturing_method / register_method!

# --- Metrics & scoring ---
export score_suite_absolute, apply_scores!, collect_hard_gate_failures
export sample_solution_1d, shock_thickness_sp, excess_dissipation, smooth_region_mask
export observed_orders, order_pass, solution_extrema, overshoot_metric
export SodProblem, sod_exact, sod_exact_conserved

# --- Verification cases (1D) — public case runners ---
export run_advection_smooth_order, run_advection_conservation
export run_burgers_conservation, run_burgers_oscillation
export burgers_square_ic
export run_euler_smooth_order, run_euler_conservation
export run_bc_transmissive_test, run_bc_dirichlet_test
export euler_density_wave_conserved
export run_persson_vs_null_burgers
export run_sod, run_shu_osher, run_m5_quant_suite, shu_osher_ic, sod_ic
# Suite orchestrators (run_m1_… / run_p31_…) are internal: FRForge.run_m1_advection_suite etc.

# --- Verification cases (2D / curved / optional) ---
export run_advection2d_smooth_order, run_euler2d_smooth_order, run_euler2d_discontinuous
export run_freestream_preservation_2d, run_advection2d_curved_order
export run_euler2d_curved_discontinuous
export riemann2d_ic, RIEMANN2D_CFG6, run_euler2d_riemann
export isentropic_vortex_primitives, run_isentropic_vortex_order
export double_mach_states, double_mach_ic, run_double_mach_reflection
export make_ffs_solid_mask, run_forward_facing_step
# run_m8/p31/p32/p33* suites, riemann2d_cfg3_ic / RIEMANN2D_CFG3: internal / research

# --- Invent / score / log / robustness / snapshot / confirm ---
export invent_method, score_reports, run_method_report
export classify_candidate, print_candidate_summary, score_suite_relative
export tradeoff_ok, DEFAULT_SCORE_MARGIN
export FROZEN_INVENT_SCHEME, default_experiment_log_path, default_experiment_log_yaml_path
export make_entry_id, entry_from_invent, format_entry_markdown
export append_experiment_entry!, invent_append_log!, list_experiment_entry_ids
export narrative_required, package_root
# NARRATIVE_PLACEHOLDER is internal; tests may use FRForge.NARRATIVE_PLACEHOLDER
export parse_experiment_log, parse_experiment_log_text, get_experiment_entry
export log_summary, log_frontier, log_pareto, log_lessons
export format_log_summary_text, format_log_frontier_text, format_log_pareto_text
export format_log_lessons_text, format_log_entry_text
export NEAR_MISS_COMPOSITE_MARGIN
export resolve_method_sources, freeze_snapshot, verify_snapshot, snapshot_tables
export SNAPSHOT_SCHEMA_VERSION
# METHOD_SOURCE_MAP is internal; use resolve_method_sources
export robustness_cells, scheme_slug, evaluate_robustness_cell, run_robustness_matrix
export assess_publication_grade, cell_ok, order_preserved, append_robustness_log_entry!
export ConfirmMeshSpec, CONFIRM_MESH_PRESETS, get_confirm_preset
export mesh_summary_dict, mesh_note_string
export run_confirm_suite, run_confirm_report, classify_confirm, confirm_method
export print_confirm_summary, entry_from_confirm, method_has_confirm_pass

# --- High-order VTK (docs / local visualization) ---
export write_vtu_high_order, write_vtu_high_order_with_capturing
export compute_capturing_diagnostics_2d
export vtk_lagrange_line_nodes, vtk_lagrange_quad_nodes
export vtk_point_counts_1d, vtk_point_counts_2d
export gl_to_equi_interp, parse_vtu_basic, VTK_LAGRANGE_LINE, VTK_LAGRANGE_QUAD
# Advanced AV / mesh helpers (not invent path; qualify as FRForge.* if needed):
# dissipation_operator, sensor_field_2d, element_viscosities_2d,
# element_max_wavespeed_2d, wavy_physical, wavy_partials, METHOD_REGISTRY,
# METHOD_SOURCE_MAP, NARRATIVE_PLACEHOLDER

end # module
