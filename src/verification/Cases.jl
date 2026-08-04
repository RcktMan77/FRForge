# Verification cases: 1D advection, Burgers, Euler, capturing, quant (Sod / Shu–Osher).

"""
    _order_study_dt(p, Ne_fine, a; cfl, t_final)

Fixed Δt for spatial order studies: based on the **finest** mesh and scaled
so SSP-RK3 temporal error does not dominate formal spatial order p+1.
"""
function _order_study_dt(p::Int, Ne_fine::Int, a::Float64; cfl::Float64 = 0.2)
    h = 1.0 / Ne_fine
    λ = abs(a)
    # SSP-RK3 is O(dt³); keep temporal error below spatial O(h^{p+1}) on the
    # finest mesh by shrinking CFL with p (and an extra safety factor).
    cfl_eff = cfl / (p + 1)
    return cfl_eff * h / ((2p + 1) * max(λ, eps(Float64)))
end

"""
    run_advection_smooth_order(; p, n_elements_list, a, t_final, cfl, order_tol)

Smooth sine-wave advection order study on periodic domain [0,1].
Uses a **fixed** Δt from the finest mesh so observed order measures space.
Returns a case dict suitable for JSON reports.
"""
function run_advection_smooth_order(;
    p::Int = 3,
    n_elements_list::AbstractVector{Int} = [8, 16, 32, 64],
    a::Float64 = 1.0,
    t_final::Float64 = 1.0,  # one period on [0,1]
    cfl::Float64 = 0.2,
    order_tol::Float64 = DEFAULT_ORDER_TOLERANCE,
    method::AbstractCapturingMethod = NullCapturing(),
)
    t0 = time()
    eq = LinearAdvection1D(a)
    ops = build_operators(p)
    formal = p + 1
    Ne_fine = maximum(n_elements_list)
    dt_fixed = _order_study_dt(p, Ne_fine, a; cfl = cfl)

    mesh_sizes = Float64[]
    l2_errors = Float64[]
    diverged = false
    nan_detected = false

    uexact = x -> sin(2π * (x - a * t_final))

    for Ne in n_elements_list
        mesh = Mesh1D(0.0, 1.0, Ne; left_bc = PeriodicBC(), right_bc = PeriodicBC())
        state = allocate_state(mesh, ops, Val(1))
        set_initial_condition!(state, x -> sin(2π * x))
        result = ssp_rk3!(state, eq, method, t_final; dt = dt_fixed)
        if result.status != :ok
            diverged = true
            nan_detected = true
            push!(mesh_sizes, 1.0 / Ne)
            push!(l2_errors, NaN)
            continue
        end
        push!(mesh_sizes, 1.0 / Ne)
        push!(l2_errors, l2_error(state, uexact, 1))
        if has_nonfinite(state.u)
            nan_detected = true
            diverged = true
        end
    end

    obs = observed_orders(mesh_sizes, l2_errors)
    opass = !diverged && order_pass(obs; formal_order = formal, tol = order_tol)

    # Conservation on finest mesh (nonzero-mean IC for a well-defined relative residual)
    Ne = n_elements_list[end]
    mesh = Mesh1D(0.0, 1.0, Ne; left_bc = PeriodicBC(), right_bc = PeriodicBC())
    state = allocate_state(mesh, ops, Val(1))
    set_initial_condition!(state, x -> 1.0 + 0.5 * sin(2π * x))
    M0 = discrete_mass(state, 1)
    result = ssp_rk3!(state, eq, method, t_final; cfl = cfl)
    MT = discrete_mass(state, 1)
    cres = conservation_residual_relative(M0, MT)
    abs_res = conservation_residual_absolute(M0, MT)
    cpass = result.status == :ok && (cres <= 1e-12 || abs_res <= 1e-12)

    case_pass = opass && cpass && !diverged && !nan_detected

    return case_report_dict(;
    name="advection_smooth_order_p$(p)",
    case_type="smooth_order",
    equation="linear_advection",
    p=p,
    capturing_method="null",
    pass=case_pass,
    diverged=diverged,
    nan_detected=nan_detected,
    conservation_residual=cres,
    conservation_pass=cpass,
    conservation_metric="periodic_mass_change",
    positivity_ok=true,
    wall_time_sec=time() - t0,
    metrics=Dict{String,Any}(
            "a" => a,
            "t_final" => t_final,
            "cfl" => cfl,
            "dt_fixed" => dt_fixed,
            "mass_abs_change" => abs_res,
            "n_elements_list" => collect(n_elements_list),
        ),
    extra=Dict{String,Any}(
        "mesh_sizes" => mesh_sizes,
        "l2_errors" => l2_errors,
        "observed_orders" => obs,
        "formal_order" => formal,
        "order_pass" => opass,
        "order_tolerance" => order_tol,
    ),
)
end

"""
    run_advection_conservation(; p, n_elements, t_final, cfl)

Dedicated periodic conservation test for linear advection.
Uses IC with nonzero mean so relative residual is well-defined.
"""
function run_advection_conservation(;
    p::Int = 3,
    n_elements::Int = 16,
    a::Float64 = 1.0,
    t_final::Float64 = 1.0,
    cfl::Float64 = 0.2,
    method::AbstractCapturingMethod = NullCapturing(),
)
    t0 = time()
    eq = LinearAdvection1D(a)
    ops = build_operators(p)
    mesh = Mesh1D(0.0, 1.0, n_elements; left_bc = PeriodicBC(), right_bc = PeriodicBC())
    state = allocate_state(mesh, ops, Val(1))
    # Nonzero mean → mass O(1); relative residual is meaningful
    set_initial_condition!(state, x -> 1.0 + 0.5 * sin(2π * x))
    M0 = discrete_mass(state, 1)
    result = ssp_rk3!(state, eq, method, t_final; cfl = cfl)
    MT = discrete_mass(state, 1)
    cres = conservation_residual_relative(M0, MT)
    abs_res = conservation_residual_absolute(M0, MT)
    # Target near machine precision (design: ~1e-12 – 1e-14 relative)
    cpass = result.status == :ok && (cres <= 1e-12 || abs_res <= 1e-13)
    diverged = result.status != :ok
    nan_detected = diverged || has_nonfinite(state.u)

    return case_report_dict(;
    name="advection_conservation_p$(p)",
    case_type="other",
    equation="linear_advection",
    p=p,
    capturing_method="null",
    pass=cpass && !diverged,
    diverged=diverged,
    nan_detected=nan_detected,
    conservation_residual=cres,
    conservation_pass=cpass,
    conservation_metric="periodic_mass_change",
    positivity_ok=true,
    wall_time_sec=time() - t0,
    metrics=Dict{String,Any}(
            "mass_initial" => M0,
            "mass_final" => MT,
            "mass_abs_change" => abs_res,
            "n_elements" => n_elements,
            "n_steps" => result.n_steps,
        ),
)
end

"""
    run_m1_advection_suite() -> (cases, overall_pass, hard_gate_failures)

Run 1D advection verification: conservation + order for p=2,3,4.
"""
function run_m1_advection_suite()
    cases = Any[]
    hard_fails = String[]

    for p in (2, 3, 4)
        # Enough refinements for robust slopes. For high p, stop before the
        # temporal/roundoff floor on ultra-fine meshes (still ≥3 grids).
        nlist = p >= 4 ? [4, 8, 16] : [8, 16, 32, 64]
        c = run_advection_smooth_order(; p = p, n_elements_list = nlist)
        push!(cases, c)
        if !c["order_pass"]
            push!(hard_fails, "order_pass failed for p=$p (observed=$(c["observed_orders"]))")
        end
        if !c["conservation_pass"]
            push!(hard_fails, "conservation_pass failed for order case p=$p")
        end

        cc = run_advection_conservation(; p = p, n_elements = 32)
        push!(cases, cc)
        if !cc["conservation_pass"]
            push!(
                hard_fails,
                "conservation residual too large for p=$p: $(cc["conservation_residual"])",
            )
        end
    end

    overall = all(c -> c["pass"] === true, cases) && isempty(hard_fails)
    return cases, overall, hard_fails
end

# ---------------------------------------------------------------------------
# --- 1D inviscid Burgers ---
# ---------------------------------------------------------------------------

"""Periodic square-wave IC for Burgers oscillation demo on [0,1].

`x_disc` should not sit exactly on an element face when using modal sensors
(otherwise each element is piecewise constant and the Persson indicator is zero).
Default 1/3 places the jump inside an element for typical even meshes.
"""
function burgers_square_ic(x; u_left = 1.0, u_right = 0.0, x_disc = 1 / 3)
    return x < x_disc ? u_left : u_right
end
"""
    run_burgers_conservation(; p, n_elements, t_final, cfl)

Periodic Burgers with discontinuous IC; check mass conservation.
"""
function run_burgers_conservation(;
    p::Int = 3,
    n_elements::Int = 32,
    t_final::Float64 = 0.15,
    cfl::Float64 = 0.2,
    method::AbstractCapturingMethod = NullCapturing(),
)
    t0 = time()
    eq = Burgers1D()
    ops = build_operators(p)
    mesh = Mesh1D(0.0, 1.0, n_elements; left_bc = PeriodicBC(), right_bc = PeriodicBC())
    state = allocate_state(mesh, ops, Val(1))
    set_initial_condition!(state, x -> burgers_square_ic(x))
    M0 = discrete_mass(state, 1)
    result = ssp_rk3!(state, eq, method, t_final; cfl = cfl)
    MT = discrete_mass(state, 1)
    cres = conservation_residual_relative(M0, MT)
    abs_res = conservation_residual_absolute(M0, MT)
    cpass = result.status == :ok && (cres <= 1e-10 || abs_res <= 1e-12)
    diverged = result.status != :ok
    nan_detected = diverged || has_nonfinite(state.u)

    return case_report_dict(;
    name="burgers_conservation_p$(p)",
    case_type="other",
    equation="burgers",
    p=p,
    capturing_method="null",
    pass=cpass && !diverged,
    diverged=diverged,
    nan_detected=nan_detected,
    conservation_residual=cres,
    conservation_pass=cpass,
    conservation_metric="periodic_mass_change",
    positivity_ok=true,
    wall_time_sec=time() - t0,
    metrics=Dict{String,Any}(
            "mass_initial" => M0,
            "mass_final" => MT,
            "mass_abs_change" => abs_res,
            "n_elements" => n_elements,
            "n_steps" => result.n_steps,
            "t_final" => t_final,
        ),
)
end

"""
    run_burgers_oscillation(; p, n_elements, t_final, cfl, min_overshoot)

Discontinuous Burgers with pure high-order FR (NullCapturing).

Success means the run completes, conserves, and **exhibits** Gibbs-type
oscillations (overshoot above IC max or undershoot below IC min). That documents
the failure mode that shock-capturing methods later aim to fix.
"""
function run_burgers_oscillation(;
    p::Int = 3,
    n_elements::Int = 32,
    t_final::Float64 = 0.15,
    cfl::Float64 = 0.2,
    min_overshoot::Float64 = 0.02,  # require at least 2% of jump height
    method::AbstractCapturingMethod = NullCapturing(),
    u_left::Float64 = 1.0,
    u_right::Float64 = 0.0,
)
    t0 = time()
    eq = Burgers1D()
    ops = build_operators(p)
    mesh = Mesh1D(0.0, 1.0, n_elements; left_bc = PeriodicBC(), right_bc = PeriodicBC())
    state = allocate_state(mesh, ops, Val(1))
    set_initial_condition!(state, x -> burgers_square_ic(x; u_left = u_left, u_right = u_right))

    u0_min, u0_max = solution_extrema(state, 1)
    M0 = discrete_mass(state, 1)

    result = ssp_rk3!(state, eq, method, t_final; cfl = cfl)
    MT = discrete_mass(state, 1)
    u_min, u_max = solution_extrema(state, 1)
    η_over, η_under, η = overshoot_metric(u_min, u_max, u0_min, u0_max)

    cres = conservation_residual_relative(M0, MT)
    abs_res = conservation_residual_absolute(M0, MT)
    cpass = result.status == :ok && (cres <= 1e-10 || abs_res <= 1e-12)
    diverged = result.status != :ok
    nan_detected = diverged || has_nonfinite(state.u)

    # Success: oscillations present (expected HO failure mode without capturing)
    oscillations_present = η >= min_overshoot
    # Case "pass" means demo succeeded: ran, conserved, and showed oscillations
    case_pass = cpass && !diverged && !nan_detected && oscillations_present

    return case_report_dict(;
    name="burgers_oscillation_p$(p)",
    case_type="discontinuous",
    equation="burgers",
    p=p,
    capturing_method="null",
    pass=case_pass,
    diverged=diverged,
    nan_detected=nan_detected,
    conservation_residual=cres,
    conservation_pass=cpass,
    conservation_metric="periodic_mass_change",
    positivity_ok=true,
    wall_time_sec=time() - t0,
    n_elements=n_elements,
    t_final=t_final,
    excess_dissipation=nothing,
    shock_thickness=nothing,
    shock_thickness_unit="sp_spacings",
    overshoot=η,
    metrics=Dict{String,Any}(
            "mass_initial" => M0,
            "mass_final" => MT,
            "mass_abs_change" => abs_res,
            "u_min" => u_min,
            "u_max" => u_max,
            "u0_min" => u0_min,
            "u0_max" => u0_max,
            "overshoot_high" => η_over,
            "overshoot_low" => η_under,
            "overshoot" => η,
            "min_overshoot_required" => min_overshoot,
            "oscillations_present" => oscillations_present,
            "n_steps" => result.n_steps,
            "note" => "Documents HO oscillatory failure of NullCapturing; oscillations are expected",
        ),
)
end

"""
    run_m2_burgers_suite() -> (cases, overall_pass, hard_gate_failures)

Burgers conservation + oscillatory HO failure for p=2,3,4.
"""
function run_m2_burgers_suite()
    cases = Any[]
    hard_fails = String[]

    for p in (2, 3, 4)
        cc = run_burgers_conservation(; p = p, n_elements = 32, t_final = 0.15)
        push!(cases, cc)
        if !cc["conservation_pass"]
            push!(
                hard_fails,
                "burgers conservation failed for p=$p: res=$(cc["conservation_residual"])",
            )
        end

        co = run_burgers_oscillation(; p = p, n_elements = 32, t_final = 0.15)
        push!(cases, co)
        if co["diverged"] || co["nan_detected"]
            push!(hard_fails, "burgers oscillation run diverged for p=$p")
        elseif !co["metrics"]["oscillations_present"]
            push!(
                hard_fails,
                "burgers p=$p did not show expected HO overshoot (η=$(co["overshoot"]))",
            )
        elseif !co["conservation_pass"]
            push!(hard_fails, "burgers oscillation conservation failed for p=$p")
        end
    end

    overall = all(c -> c["pass"] === true, cases) && isempty(hard_fails)
    return cases, overall, hard_fails
end

# ---------------------------------------------------------------------------
# --- 1D Euler + BCs + smooth order ---
# ---------------------------------------------------------------------------

"""
Exact density-wave solution (periodic):
  ρ = 1 + A sin(2π(x - u0 t)),  u = u0,  p = p0  (constant)
Conserved: (ρ, ρu, E) with E = p/(γ-1) + ½ ρ u².
"""
function euler_density_wave_conserved(
    eq::Euler1D{T},
    x,
    t;
    A = T(0.2),
    u0 = T(1),
    p0 = T(1),
) where {T}
    ρ = one(T) + T(A) * sin(2π * (x - u0 * t))
    return primitives_to_conserved(eq, ρ, u0, p0)
end

function _euler_order_study_dt(p::Int, Ne_fine::Int, λ_max::Float64; cfl::Float64 = 0.2)
    h = 1.0 / Ne_fine
    cfl_eff = cfl / (p + 1)
    return cfl_eff * h / ((2p + 1) * max(λ_max, eps(Float64)))
end

"""Positivity over a full SolutionState for Euler."""
function positivity_ok(eq::Euler1D, state::SolutionState; kwargs...)
    return positivity_ok(eq, state.u; kwargs...)
end

"""
    run_euler_smooth_order(; p, n_elements_list, t_final, cfl, order_tol)

Smooth density-wave order study on periodic [0,1] with NullCapturing.
Formal order target p+1 (density L2).
"""
function run_euler_smooth_order(;
    p::Int = 3,
    n_elements_list::AbstractVector{Int} = [8, 16, 32],
    t_final::Float64 = 1.0,
    cfl::Float64 = 0.2,
    order_tol::Float64 = DEFAULT_ORDER_TOLERANCE,
    γ::Float64 = 1.4,
    method::AbstractCapturingMethod = NullCapturing(),
    scheme::SchemeConfig = DEFAULT_SCHEME,
)
    t0 = time()
    eq = Euler1D(γ)
    ops = build_operators(p; points = scheme.points)
    formal = p + 1
    Ne_fine = maximum(n_elements_list)
    λ_max = 2.5
    dt_fixed = _euler_order_study_dt(p, Ne_fine, λ_max; cfl = cfl)

    mesh_sizes = Float64[]
    l2_errors = Float64[]
    diverged = false
    nan_detected = false
    pos_ok = true

    uexact = x -> euler_density_wave_conserved(eq, x, t_final)

    for Ne in n_elements_list
        mesh = Mesh1D(0.0, 1.0, Ne; left_bc = PeriodicBC(), right_bc = PeriodicBC())
        state = allocate_state(mesh, ops, Val(3); scheme = scheme)
        set_initial_condition!(state, x -> euler_density_wave_conserved(eq, x, 0.0))
        result = integrate!(state, eq, method, t_final; dt = dt_fixed)
        if result.status != :ok
            diverged = true
            nan_detected = true
            push!(mesh_sizes, 1.0 / Ne)
            push!(l2_errors, NaN)
            continue
        end
        push!(mesh_sizes, 1.0 / Ne)
        push!(l2_errors, l2_error(state, uexact, 1))
        pos_ok = pos_ok && positivity_ok(eq, state)
        if has_nonfinite(state.u)
            nan_detected = true
            diverged = true
        end
    end

    obs = observed_orders(mesh_sizes, l2_errors)
    opass = !diverged && order_pass(obs; formal_order = formal, tol = order_tol)

    Ne = n_elements_list[end]
    mesh = Mesh1D(0.0, 1.0, Ne; left_bc = PeriodicBC(), right_bc = PeriodicBC())
    state = allocate_state(mesh, ops, Val(3); scheme = scheme)
    set_initial_condition!(state, x -> euler_density_wave_conserved(eq, x, 0.0))
    M0 = [discrete_mass(state, c) for c in 1:3]
    result = integrate!(state, eq, method, t_final; cfl = cfl)
    MT = [discrete_mass(state, c) for c in 1:3]
    cres_comp = [conservation_residual_relative(M0[c], MT[c]) for c in 1:3]
    abs_comp = [conservation_residual_absolute(M0[c], MT[c]) for c in 1:3]
    cres = maximum(cres_comp)
    cpass =
        result.status == :ok &&
        all(cres_comp[c] <= 1e-10 || abs_comp[c] <= 1e-11 for c in 1:3)
    pos_ok = pos_ok && positivity_ok(eq, state)

    case_pass = opass && cpass && !diverged && !nan_detected && pos_ok

    return case_report_dict(;
    name="euler_smooth_order_p$(p)",
    case_type="smooth_order",
    equation="euler1d",
    p=p,
    capturing_method="null",
    pass=case_pass,
    diverged=diverged,
    nan_detected=nan_detected,
    conservation_residual=cres,
    conservation_pass=cpass,
    conservation_metric="periodic_mass_change",
    positivity_ok=pos_ok,
    wall_time_sec=time() - t0,
    metrics=Dict{String,Any}(
            "γ" => γ,
            "t_final" => t_final,
            "cfl" => cfl,
            "dt_fixed" => dt_fixed,
            "conservation_residual_components" => cres_comp,
            "mass_abs_change_components" => abs_comp,
            "n_elements_list" => collect(n_elements_list),
        ),
    extra=Dict{String,Any}(
        "mesh_sizes" => mesh_sizes,
        "l2_errors" => l2_errors,
        "observed_orders" => obs,
        "formal_order" => formal,
        "order_pass" => opass,
        "order_tolerance" => order_tol,
    ),
)
end

"""
    run_euler_conservation(; p, n_elements, t_final, cfl)

Periodic Euler density wave: mass/momentum/energy conservation.
"""
function run_euler_conservation(;
    p::Int = 3,
    n_elements::Int = 32,
    t_final::Float64 = 1.0,
    cfl::Float64 = 0.2,
    γ::Float64 = 1.4,
    method::AbstractCapturingMethod = NullCapturing(),
)
    t0 = time()
    eq = Euler1D(γ)
    ops = build_operators(p)
    mesh = Mesh1D(0.0, 1.0, n_elements; left_bc = PeriodicBC(), right_bc = PeriodicBC())
    state = allocate_state(mesh, ops, Val(3))
    set_initial_condition!(state, x -> euler_density_wave_conserved(eq, x, 0.0))
    M0 = [discrete_mass(state, c) for c in 1:3]
    result = ssp_rk3!(state, eq, method, t_final; cfl = cfl)
    MT = [discrete_mass(state, c) for c in 1:3]
    cres_comp = [conservation_residual_relative(M0[c], MT[c]) for c in 1:3]
    abs_comp = [conservation_residual_absolute(M0[c], MT[c]) for c in 1:3]
    cres = maximum(cres_comp)
    cpass =
        result.status == :ok &&
        all(cres_comp[c] <= 1e-11 || abs_comp[c] <= 1e-12 for c in 1:3)
    diverged = result.status != :ok
    nan_detected = diverged || has_nonfinite(state.u)
    pos = positivity_ok(eq, state)

    return case_report_dict(;
    name="euler_conservation_p$(p)",
    case_type="other",
    equation="euler1d",
    p=p,
    capturing_method="null",
    pass=cpass && !diverged && pos,
    diverged=diverged,
    nan_detected=nan_detected,
    conservation_residual=cres,
    conservation_pass=cpass,
    conservation_metric="periodic_mass_change",
    positivity_ok=pos,
    wall_time_sec=time() - t0,
    metrics=Dict{String,Any}(
            "mass_initial" => M0,
            "mass_final" => MT,
            "mass_abs_change_components" => abs_comp,
            "conservation_residual_components" => cres_comp,
            "n_elements" => n_elements,
            "n_steps" => result.n_steps,
        ),
)
end

"""
    run_bc_transmissive_test(; p, n_elements, t_final)

Uniform free-stream with TransmissiveBC both ends — state should stay uniform.
"""
function run_bc_transmissive_test(;
    p::Int = 3,
    n_elements::Int = 16,
    t_final::Float64 = 0.2,
    cfl::Float64 = 0.2,
    γ::Float64 = 1.4,
)
    t0 = time()
    eq = Euler1D(γ)
    ops = build_operators(p)
    mesh = Mesh1D(
        0.0,
        1.0,
        n_elements;
        left_bc = TransmissiveBC(),
        right_bc = TransmissiveBC(),
    )
    state = allocate_state(mesh, ops, Val(3))
    U0 = primitives_to_conserved(eq, 1.0, 0.5, 1.0)
    set_initial_condition!(state, _ -> U0)
    result = ssp_rk3!(state, eq, NullCapturing(), t_final; cfl = cfl)
    dev = 0.0
    for c in 1:3
        for e in 1:n_elements, j in 1:size(state.u, 1)
            dev = max(dev, abs(state.u[j, e, c] - U0[c]))
        end
    end
    ok = result.status == :ok && dev < 1e-10 && positivity_ok(eq, state)

    return case_report_dict(;
    name="bc_transmissive_freestream_p$(p)",
    case_type="other",
    equation="euler1d",
    p=p,
    capturing_method="null",
    pass=ok,
    diverged=result.status != :ok,
    nan_detected=has_nonfinite(state.u),
    conservation_residual=dev,
    conservation_pass=dev < 1e-10,
    conservation_metric="none",
    positivity_ok=positivity_ok(eq, state),
    wall_time_sec=time() - t0,
    metrics=Dict{String,Any}(
            "max_deviation_from_uniform" => dev,
            "n_steps" => result.n_steps,
            "bc" => "transmissive",
        ),
)
end

"""
    run_bc_dirichlet_test(; p, n_elements, t_final)

Dirichlet freestream both ends with uniform IC — should remain uniform.
"""
function run_bc_dirichlet_test(;
    p::Int = 3,
    n_elements::Int = 16,
    t_final::Float64 = 0.2,
    cfl::Float64 = 0.2,
    γ::Float64 = 1.4,
)
    t0 = time()
    eq = Euler1D(γ)
    ops = build_operators(p)
    U0 = primitives_to_conserved(eq, 1.0, 0.3, 1.0)
    mesh = Mesh1D(
        0.0,
        1.0,
        n_elements;
        left_bc = DirichletBC(t -> U0),
        right_bc = DirichletBC(t -> U0),
    )
    state = allocate_state(mesh, ops, Val(3))
    set_initial_condition!(state, _ -> U0)
    result = ssp_rk3!(state, eq, NullCapturing(), t_final; cfl = cfl)
    dev = 0.0
    for c in 1:3
        for e in 1:n_elements, j in 1:size(state.u, 1)
            dev = max(dev, abs(state.u[j, e, c] - U0[c]))
        end
    end
    ok = result.status == :ok && dev < 1e-10 && positivity_ok(eq, state)

    return case_report_dict(;
    name="bc_dirichlet_freestream_p$(p)",
    case_type="other",
    equation="euler1d",
    p=p,
    capturing_method="null",
    pass=ok,
    diverged=result.status != :ok,
    nan_detected=has_nonfinite(state.u),
    conservation_residual=dev,
    conservation_pass=dev < 1e-10,
    conservation_metric="none",
    positivity_ok=positivity_ok(eq, state),
    wall_time_sec=time() - t0,
    metrics=Dict{String,Any}(
            "max_deviation_from_uniform" => dev,
            "n_steps" => result.n_steps,
            "bc" => "dirichlet",
        ),
)
end

"""
    run_m3_euler_suite() -> (cases, overall_pass, hard_gate_failures)

Euler smooth order (p=2,3,4), conservation, and BC path tests.
"""
function run_m3_euler_suite()
    cases = Any[]
    hard_fails = String[]

    for p in (2, 3, 4)
        # Pre-asymptotic rates on very coarse grids; use meshes in the asymptotic regime
        nlist = p == 2 ? [16, 32, 64] : (p == 3 ? [8, 16, 32] : [16, 32, 64])
        c = run_euler_smooth_order(; p = p, n_elements_list = nlist, cfl = 0.1)
        push!(cases, c)
        if !c["order_pass"]
            push!(hard_fails, "euler order failed p=$p observed=$(c["observed_orders"])")
        end
        if !c["conservation_pass"]
            push!(hard_fails, "euler order-case conservation failed p=$p")
        end
        if !c["positivity_ok"]
            push!(hard_fails, "euler positivity failed p=$p")
        end

        cc = run_euler_conservation(; p = p, n_elements = 32)
        push!(cases, cc)
        if !cc["conservation_pass"]
            push!(hard_fails, "euler conservation failed p=$p res=$(cc["conservation_residual"])")
        end
    end

    bt = run_bc_transmissive_test()
    push!(cases, bt)
    if !bt["pass"]
        push!(hard_fails, "transmissive BC freestream test failed")
    end

    bd = run_bc_dirichlet_test()
    push!(cases, bd)
    if !bd["pass"]
        push!(hard_fails, "dirichlet BC freestream test failed")
    end

    overall = all(c -> c["pass"] === true, cases) && isempty(hard_fails)
    return cases, overall, hard_fails
end

# ---------------------------------------------------------------------------
# --- Capturing interface + Persson AV baseline ---
# ---------------------------------------------------------------------------

"""
    run_persson_vs_null_burgers(; p, n_elements, t_final)

Compare NullCapturing vs PerssonAV on discontinuous Burgers.
Success: both run; Persson reduces overshoot vs null; conservation holds for both.
"""
function run_persson_vs_null_burgers(;
    p::Int = 3,
    n_elements::Int = 32,
    t_final::Float64 = 0.15,
    cfl::Float64 = 0.2,
)
    t0 = time()
    c_null = run_burgers_oscillation(;
        p = p,
        n_elements = n_elements,
        t_final = t_final,
        cfl = cfl,
        method = NullCapturing(),
    )
    # Rename for clarity in report
    c_null = deepcopy(c_null)
    c_null["name"] = "burgers_null_p$(p)"
    c_null["capturing_method"] = "null"
    c_null["method_params"] = method_params(NullCapturing())

    eq = Burgers1D()
    ops = build_operators(p)
    mesh = Mesh1D(0.0, 1.0, n_elements; left_bc = PeriodicBC(), right_bc = PeriodicBC())
    state = allocate_state(mesh, ops, Val(1))
    set_initial_condition!(state, x -> burgers_square_ic(x))
    u0_min, u0_max = solution_extrema(state, 1)
    M0 = discrete_mass(state, 1)
    method = PerssonAVMethod()
    # AV stiffens the residual; use a slightly tighter CFL than pure hyperbolic
    result = ssp_rk3!(state, eq, method, t_final; cfl = min(cfl, 0.1))
    MT = discrete_mass(state, 1)
    u_min, u_max = solution_extrema(state, 1)
    _, _, η = overshoot_metric(u_min, u_max, u0_min, u0_max)
    cres = conservation_residual_relative(M0, MT)
    abs_res = conservation_residual_absolute(M0, MT)
    cpass = result.status == :ok && (cres <= 1e-9 || abs_res <= 1e-11)
    diverged = result.status != :ok
    nan_detected = diverged || has_nonfinite(state.u)
    η_null = c_null["overshoot"]
    reduced = η < η_null  # Persson should damp oscillations relative to null
    # Soft requirement: either reduced overshoot OR sensor activated (σ mean > 0)
    # Hard: run stable + conserve
    case_pass = cpass && !diverged && !nan_detected

    c_pers = case_report_dict(;
    name="burgers_persson_av_p$(p)",
    case_type="discontinuous",
    equation="burgers",
    p=p,
    capturing_method="persson_av",
    pass=case_pass,
    diverged=diverged,
    nan_detected=nan_detected,
    conservation_residual=cres,
    conservation_pass=cpass,
    conservation_metric="periodic_mass_change",
    positivity_ok=true,
    wall_time_sec=time() - t0,
    n_elements=n_elements,
    t_final=t_final,
    excess_dissipation=nothing,
    shock_thickness=nothing,
    shock_thickness_unit="sp_spacings",
    overshoot=η,
    metrics=Dict{String,Any}(
            "overshoot" => η,
            "overshoot_null" => η_null,
            "overshoot_reduced_vs_null" => reduced,
            "mass_abs_change" => abs_res,
            "n_steps" => result.n_steps,
            "av_form" => method.dissip.av_form,
        ),
    extra=Dict{String,Any}(
        "method_params" => method_params(method),
    ),
)

    # Comparison summary case
    both_ok = c_null["pass"] && c_pers["pass"] && reduced
    c_cmp = case_report_dict(;
    name="burgers_persson_vs_null_p$(p)",
    case_type="other",
    equation="burgers",
    p=p,
    capturing_method="persson_av",
    pass=both_ok,
    diverged=false,
    nan_detected=false,
    conservation_residual=max(c_null["conservation_residual"], c_pers["conservation_residual"]),
    conservation_pass=c_null["conservation_pass"] && c_pers["conservation_pass"],
    conservation_metric="periodic_mass_change",
    positivity_ok=true,
    wall_time_sec=time() - t0,
    metrics=Dict{String,Any}(
            "overshoot_null" => η_null,
            "overshoot_persson" => η,
            "overshoot_reduced" => reduced,
            "note" => "Persson AV baseline: should reduce HO overshoot vs NullCapturing",
        ),
)

    return c_null, c_pers, c_cmp
end

"""
    run_m4_capturing_suite() -> (cases, overall_pass, hard_gate_failures)
"""
function run_m4_capturing_suite()
    cases = Any[]
    hard_fails = String[]

    # Viscous operator mass residual (constants) on a sample mesh
    ops = build_operators(3)
    mesh = Mesh1D(0.0, 1.0, 8)
    mass_res = viscous_mass_residual_scale(ops, mesh)
    c_mass = case_report_dict(;
    name="av_conservative_br0_mass_check",
    case_type="other",
    equation="operator",
    p=3,
    capturing_method="persson_av",
    pass=mass_res < 1e-12,
    diverged=false,
    nan_detected=false,
    conservation_residual=mass_res,
    conservation_pass=mass_res < 1e-12,
    conservation_metric="periodic_mass_change",
    positivity_ok=true,
    wall_time_sec=0.0,
    metrics=Dict{String,Any}(
            "weighted_mass_of_AV_on_constant" => mass_res,
            "av_form" => "conservative_br0",
        ),
)
    push!(cases, c_mass)
    if !c_mass["pass"]
        push!(
            hard_fails,
            "conservative_br0 AV does not annihilate constants in mass sense: $mass_res",
        )
    end

    for p in (2, 3, 4)
        c_null, c_pers, c_cmp = run_persson_vs_null_burgers(; p = p)
        push!(cases, c_null, c_pers, c_cmp)
        if !c_pers["conservation_pass"]
            push!(hard_fails, "PerssonAV conservation failed p=$p")
        end
        if c_pers["diverged"]
            push!(hard_fails, "PerssonAV diverged p=$p")
        end
        if !c_cmp["metrics"]["overshoot_reduced"]
            push!(
                hard_fails,
                "PerssonAV did not reduce overshoot vs null for p=$p (null=$(c_cmp["metrics"]["overshoot_null"]), pers=$(c_cmp["metrics"]["overshoot_persson"]))",
            )
        end
    end

    overall = all(c -> c["pass"] === true, cases) && isempty(hard_fails)
    return cases, overall, hard_fails
end

# ---------------------------------------------------------------------------
# --- Quantitative suite: Sod, Shu–Osher, scored summary ---
# ---------------------------------------------------------------------------

"""Shu–Osher IC on [0, 10] (classic [-5,5] shifted by +5).
Left state for x < 1; right density 1+0.2 sin(5(x-5)), u=0, p=1.
"""
function shu_osher_ic(eq::Euler1D, x; x0 = 1.0)
    if x < x0
        return primitives_to_conserved(eq, 3.857143, 2.629369, 10.33333)
    else
        xc = x - 5.0  # map [0,10] → classic [-5,5]
        ρ = 1.0 + 0.2 * sin(5 * xc)
        return primitives_to_conserved(eq, ρ, 0.0, 1.0)
    end
end

function sod_ic(eq::Euler1D, x; x0 = 0.5)
    if x < x0
        return primitives_to_conserved(eq, 1.0, 0.0, 1.0)
    else
        return primitives_to_conserved(eq, 0.125, 0.0, 0.1)
    end
end

"""
    run_sod(; p, n_elements, t_final, method, cfl)

Sod shock tube with TransmissiveBC. Metrics vs exact solution and vs NullCapturing.
"""
function run_sod(;
    p::Int = 2,
    n_elements::Int = 64,
    t_final::Float64 = 0.2,
    cfl::Float64 = 0.2,
    γ::Float64 = 1.4,
    method::AbstractCapturingMethod = PerssonAVMethod(),
    method_name::AbstractString = "persson_av",
    scheme::SchemeConfig = DEFAULT_SCHEME,
)
    t0 = time()
    eq = Euler1D(γ)
    ops = build_operators(p; points = scheme.points)
    mesh = Mesh1D(
        0.0,
        1.0,
        n_elements;
        left_bc = TransmissiveBC(),
        right_bc = TransmissiveBC(),
    )
    state = allocate_state(mesh, ops, Val(3); scheme = scheme)
    set_initial_condition!(state, x -> sod_ic(eq, x))
    u0_min, u0_max = solution_extrema(state, 1)
    M0 = discrete_mass(state, 1)

    result = integrate!(state, eq, method, t_final; cfl = min(cfl, 0.15))
    diverged = result.status != :ok
    nan_detected = diverged || has_nonfinite(state.u)
    pos = !nan_detected && positivity_ok(eq, state)
    MT = discrete_mass(state, 1)
    # Open BC: mass not conserved; report change as diagnostic only
    cres = conservation_residual_absolute(M0, MT)
    cpass = true  # Sod open domain: conservation not scored; freestream BC tests cover telescoping

    x, ρ = sample_solution_1d(state, eq; component = :density)
    umin, umax = solution_extrema(state, 1)
    _, _, η = overshoot_metric(umin, umax, 0.125, 1.0)  # Sod density bounds

    # Exact density for L1 and excess
    ρ_ex = [sod_exact(SodProblem(; γ = γ, x0 = 0.5), xi, t_final)[1] for xi in x]
    mask = smooth_region_mask(
        x,
        ρ_ex;
        grad_frac = 0.15,
        exclude_windows = [(0.45, 0.55), (0.65, 0.95)],  # contact + shock neighborhoods (t=0.2)
    )
    Dex_exact = excess_dissipation(x, ρ, ρ_ex; mask = mask)

    # NullCapturing reference on same mesh / scheme
    state_n = allocate_state(mesh, ops, Val(3); scheme = scheme)
    set_initial_condition!(state_n, x -> sod_ic(eq, x))
    rn = integrate!(state_n, eq, NullCapturing(), t_final; cfl = min(cfl, 0.15))
    Dex_null = nothing
    if rn.status == :ok
        _, ρn = sample_solution_1d(state_n, eq; component = :density)
        mask_n = smooth_region_mask(
            x,
            ρn;
            grad_frac = 0.15,
            exclude_windows = [(0.45, 0.55), (0.65, 0.95)],
        )
        Dex_null = excess_dissipation(x, ρ, ρn; mask = mask_n)
    end
    Dex = Dex_null === nothing ? Dex_exact : Dex_null

    δ = shock_thickness_sp(x, ρ)  # local jump around steepest gradient

    # L1 vs exact (all points)
    l1 = sum(abs.(ρ .- ρ_ex)) / (sum(abs.(ρ_ex)) + 1e-30)

    # Pass: ran, positive, no nan — discontinuous cases don't require order
    case_pass = !diverged && !nan_detected && pos

    return case_report_dict(;
    name="sod_p$(p)_ne$(n_elements)_$(method_name)",
    case_type="discontinuous",
    equation="euler1d",
    p=p,
    capturing_method=String(method_name),
    pass=case_pass,
    diverged=diverged,
    nan_detected=nan_detected,
    conservation_residual=cres,
    conservation_pass=cpass,
    conservation_metric="none",
    positivity_ok=pos,
    wall_time_sec=time() - t0,
    n_elements=n_elements,
    t_final=t_final,
    excess_dissipation=Dex,
    shock_thickness=δ,
    shock_thickness_unit="sp_spacings",
    overshoot=η,
    metrics=Dict{String,Any}(
            "excess_dissipation_vs_null" => Dex_null,
            "excess_dissipation_vs_exact" => Dex_exact,
            "l1_density_vs_exact" => l1,
            "mass_change" => MT - M0,
            "n_steps" => result.n_steps,
            "u_min" => umin,
            "u_max" => umax,
        ),
    extra=Dict{String,Any}(
        "l1_error_vs_reference" => l1,
        "method_params" => method_params(method),
    ),
)
end

"""
    run_shu_osher(; p, n_elements, t_final, method)

Shu–Osher problem on [0,10], TransmissiveBC, t=1.8.
Excess dissipation vs same-mesh NullCapturing; optional frozen CSV reference.
"""
function run_shu_osher(;
    p::Int = 2,
    n_elements::Int = 80,
    t_final::Float64 = 1.8,
    cfl::Float64 = 0.15,
    γ::Float64 = 1.4,
    method::AbstractCapturingMethod = PerssonAVMethod(),
    method_name::AbstractString = "persson_av",
    scheme::SchemeConfig = DEFAULT_SCHEME,
)
    t0 = time()
    eq = Euler1D(γ)
    ops = build_operators(p; points = scheme.points)
    mesh = Mesh1D(
        0.0,
        10.0,
        n_elements;
        left_bc = TransmissiveBC(),
        right_bc = TransmissiveBC(),
    )
    state = allocate_state(mesh, ops, Val(3); scheme = scheme)
    set_initial_condition!(state, x -> shu_osher_ic(eq, x))
    umin0, umax0 = solution_extrema(state, 1)
    M0 = discrete_mass(state, 1)

    result = integrate!(state, eq, method, t_final; cfl = cfl)
    diverged = result.status != :ok
    nan_detected = diverged || has_nonfinite(state.u)
    pos = !nan_detected && positivity_ok(eq, state)
    MT = discrete_mass(state, 1)
    cres = conservation_residual_absolute(M0, MT)

    x, ρ = sample_solution_1d(state, eq; component = :density)
    umin, umax = solution_extrema(state, 1)
    # Post-shock density scale ~4 for overshoot bound
    _, _, η = overshoot_metric(umin, umax, umin0, max(umax0, 4.5))

    # Null reference same mesh / scheme
    state_n = allocate_state(mesh, ops, Val(3); scheme = scheme)
    set_initial_condition!(state_n, x -> shu_osher_ic(eq, x))
    rn = integrate!(state_n, eq, NullCapturing(), t_final; cfl = cfl)
    Dex = nothing
    if rn.status == :ok
        xn, ρn = sample_solution_1d(state_n, eq; component = :density)
        mask = smooth_region_mask(
            x,
            ρn;
            grad_frac = 0.12,
            exclude_windows = [(1.5, 3.0)],  # main shock region (rough)
        )
        ρn_i = interp_1d(xn, ρn, x)
        Dex = excess_dissipation(x, ρ, ρn_i; mask = mask)
    end

    δ = shock_thickness_sp(x, ρ)

    # Optional frozen high-res reference L1
    l1_ref = nothing
    ref_path = joinpath(@__DIR__, "..", "..", "test", "data", "shu_osher_ref.csv")
    if isfile(ref_path)
        xref, ρref = load_shu_osher_reference(ref_path)
        ρref_i = interp_1d(xref, ρref, x)
        l1_ref = sum(abs.(ρ .- ρref_i)) / (sum(abs.(ρref_i)) + 1e-30)
    end

    case_pass = !diverged && !nan_detected && pos

    return case_report_dict(;
    name="shu_osher_p$(p)_ne$(n_elements)_$(method_name)",
    case_type="discontinuous",
    equation="euler1d",
    p=p,
    capturing_method=String(method_name),
    pass=case_pass,
    diverged=diverged,
    nan_detected=nan_detected,
    conservation_residual=cres,
    conservation_pass=true,
    conservation_metric="none",
    positivity_ok=pos,
    wall_time_sec=time() - t0,
    n_elements=n_elements,
    t_final=t_final,
    excess_dissipation=Dex,
    shock_thickness=δ,
    shock_thickness_unit="sp_spacings",
    overshoot=η,
    metrics=Dict{String,Any}(
            "mass_change" => MT - M0,
            "n_steps" => result.n_steps,
            "u_min" => umin,
            "u_max" => umax,
            "reference_file" => isfile(ref_path) ? "test/data/shu_osher_ref.csv" : nothing,
        ),
    extra=Dict{String,Any}(
        "l1_error_vs_reference" => l1_ref,
        "method_params" => method_params(method),
    ),
)
end

function load_shu_osher_reference(path::AbstractString)
    xs = Float64[]
    ρs = Float64[]
    open(path, "r") do io
        for line in eachline(io)
            startswith(line, "#") && continue
            isempty(strip(line)) && continue
            parts = split(line, ',')
            length(parts) >= 2 || continue
            push!(xs, parse(Float64, parts[1]))
            push!(ρs, parse(Float64, parts[2]))
        end
    end
    return xs, ρs
end

"""
    run_m5_quant_suite(; method_name="persson_av", scheme=DEFAULT_SCHEME, light=false)

Quantitative suite: Euler smooth order (p=2,3) + Sod + Shu–Osher with method.
`light=true` is the CI / robustness-matrix reduced suite (smaller meshes, short Shu–Osher).
Returns cases, overall_pass, hard_gate_failures.
"""
function run_m5_quant_suite(;
    method_name::AbstractString = "persson_av",
    scheme::SchemeConfig = DEFAULT_SCHEME,
    light::Bool = false,
)
    cases = Any[]
    method = get_capturing_method(method_name)

    if light
        # CI-light: one order study, coarse Sod, short Shu–Osher
        push!(
            cases,
            run_euler_smooth_order(;
                p = 2,
                n_elements_list = [16, 32],
                cfl = 0.1,
                scheme = scheme,
            ),
        )
        push!(
            cases,
            run_sod(;
                p = 2,
                n_elements = 32,
                t_final = 0.15,
                method = method,
                method_name = method_name,
                scheme = scheme,
            ),
        )
        push!(
            cases,
            run_shu_osher(;
                p = 1,
                n_elements = 40,
                t_final = 0.4,
                cfl = 0.1,
                method = method isa PerssonAVMethod ?
                         PerssonAVMethod(c_av = max(method.dissip.c_av, 0.3)) : method,
                method_name = method_name,
                scheme = scheme,
            ),
        )
    else
        # Smooth order (subset for runtime: p=2,3)
        for p in (2, 3)
            nlist = p == 2 ? [16, 32, 64] : [8, 16, 32]
            c = run_euler_smooth_order(; p = p, n_elements_list = nlist, cfl = 0.1, scheme = scheme)
            push!(cases, c)
        end

        push!(
            cases,
            run_sod(;
                p = 2,
                n_elements = 64,
                t_final = 0.2,
                method = method,
                method_name = method_name,
                scheme = scheme,
            ),
        )

        # Shu–Osher: p=1 is the robust HO-start setting for explicit FR+AV on this problem
        push!(
            cases,
            run_shu_osher(;
                p = 1,
                n_elements = 100,
                t_final = 1.8,
                cfl = 0.1,
                method = method isa PerssonAVMethod ?
                         PerssonAVMethod(c_av = max(method.dissip.c_av, 0.3)) : method,
                method_name = method_name,
                scheme = scheme,
            ),
        )
    end

    hard_fails = collect_hard_gate_failures(cases)
    overall = all(c -> c["pass"] === true, cases) && isempty(hard_fails)
    return cases, overall, hard_fails
end
