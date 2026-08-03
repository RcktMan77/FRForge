# Verification cases for Milestone 1 (linear advection).

"""
    _order_study_dt(p, Ne_fine, a; cfl, t_final)

Fixed Δt for spatial order studies: based on the **finest** mesh and scaled
so SSP-RK3 temporal error does not dominate formal spatial order p+1.
"""
function _order_study_dt(p::Int, Ne_fine::Int, a::Float64; cfl::Float64=0.2)
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
    p::Int=3,
    n_elements_list::AbstractVector{Int}=[8, 16, 32, 64],
    a::Float64=1.0,
    t_final::Float64=1.0,  # one period on [0,1]
    cfl::Float64=0.2,
    order_tol::Float64=DEFAULT_ORDER_TOLERANCE,
    method::AbstractCapturingMethod=NullCapturing(),
)
    t0 = time()
    eq = LinearAdvection1D(a)
    ops = build_operators(p)
    formal = p + 1
    Ne_fine = maximum(n_elements_list)
    dt_fixed = _order_study_dt(p, Ne_fine, a; cfl=cfl)

    mesh_sizes = Float64[]
    l2_errors = Float64[]
    diverged = false
    nan_detected = false

    uexact = x -> sin(2π * (x - a * t_final))

    for Ne in n_elements_list
        mesh = Mesh1D(0.0, 1.0, Ne; left_bc=PeriodicBC(), right_bc=PeriodicBC())
        state = allocate_state(mesh, ops, Val(1))
        set_initial_condition!(state, x -> sin(2π * x))
        result = ssp_rk3!(state, eq, method, t_final; dt=dt_fixed)
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
    opass = !diverged && order_pass(obs; formal_order=formal, tol=order_tol)

    # Conservation on finest mesh (nonzero-mean IC for a well-defined relative residual)
    Ne = n_elements_list[end]
    mesh = Mesh1D(0.0, 1.0, Ne; left_bc=PeriodicBC(), right_bc=PeriodicBC())
    state = allocate_state(mesh, ops, Val(1))
    set_initial_condition!(state, x -> 1.0 + 0.5 * sin(2π * x))
    M0 = discrete_mass(state, 1)
    result = ssp_rk3!(state, eq, method, t_final; cfl=cfl)
    MT = discrete_mass(state, 1)
    cres = conservation_residual_relative(M0, MT)
    abs_res = conservation_residual_absolute(M0, MT)
    cpass = result.status == :ok && (cres <= 1e-12 || abs_res <= 1e-12)

    case_pass = opass && cpass && !diverged && !nan_detected

    return Dict{String,Any}(
        "name" => "advection_smooth_order_p$(p)",
        "case_type" => "smooth_order",
        "equation" => "linear_advection",
        "p" => p,
        "capturing_method" => "null",
        "pass" => case_pass,
        "diverged" => diverged,
        "nan_detected" => nan_detected,
        "conservation_residual" => cres,
        "conservation_pass" => cpass,
        "conservation_metric" => "periodic_mass_change",
        "positivity_ok" => true,
        "wall_time_sec" => time() - t0,
        "metrics" => Dict{String,Any}(
            "a" => a,
            "t_final" => t_final,
            "cfl" => cfl,
            "dt_fixed" => dt_fixed,
            "mass_abs_change" => abs_res,
            "n_elements_list" => collect(n_elements_list),
        ),
        "mesh_sizes" => mesh_sizes,
        "l2_errors" => l2_errors,
        "observed_orders" => obs,
        "formal_order" => formal,
        "order_pass" => opass,
        "order_tolerance" => order_tol,
    )
end

"""
    run_advection_conservation(; p, n_elements, t_final, cfl)

Dedicated periodic conservation test for linear advection.
Uses IC with nonzero mean so relative residual is well-defined.
"""
function run_advection_conservation(;
    p::Int=3,
    n_elements::Int=16,
    a::Float64=1.0,
    t_final::Float64=1.0,
    cfl::Float64=0.2,
    method::AbstractCapturingMethod=NullCapturing(),
)
    t0 = time()
    eq = LinearAdvection1D(a)
    ops = build_operators(p)
    mesh = Mesh1D(0.0, 1.0, n_elements; left_bc=PeriodicBC(), right_bc=PeriodicBC())
    state = allocate_state(mesh, ops, Val(1))
    # Nonzero mean → mass O(1); relative residual is meaningful
    set_initial_condition!(state, x -> 1.0 + 0.5 * sin(2π * x))
    M0 = discrete_mass(state, 1)
    result = ssp_rk3!(state, eq, method, t_final; cfl=cfl)
    MT = discrete_mass(state, 1)
    cres = conservation_residual_relative(M0, MT)
    abs_res = conservation_residual_absolute(M0, MT)
    # Target near machine precision (design: ~1e-12 – 1e-14 relative)
    cpass = result.status == :ok && (cres <= 1e-12 || abs_res <= 1e-13)
    diverged = result.status != :ok
    nan_detected = diverged || has_nonfinite(state.u)

    return Dict{String,Any}(
        "name" => "advection_conservation_p$(p)",
        "case_type" => "other",
        "equation" => "linear_advection",
        "p" => p,
        "capturing_method" => "null",
        "pass" => cpass && !diverged,
        "diverged" => diverged,
        "nan_detected" => nan_detected,
        "conservation_residual" => cres,
        "conservation_pass" => cpass,
        "conservation_metric" => "periodic_mass_change",
        "positivity_ok" => true,
        "wall_time_sec" => time() - t0,
        "metrics" => Dict{String,Any}(
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

Run M1 verification: conservation + order for p=2,3,4.
"""
function run_m1_advection_suite()
    cases = Any[]
    hard_fails = String[]

    for p in (2, 3, 4)
        # Enough refinements for robust slopes. For high p, stop before the
        # temporal/roundoff floor on ultra-fine meshes (still ≥3 grids).
        nlist = p >= 4 ? [4, 8, 16] : [8, 16, 32, 64]
        c = run_advection_smooth_order(; p=p, n_elements_list=nlist)
        push!(cases, c)
        if !c["order_pass"]
            push!(hard_fails, "order_pass failed for p=$p (observed=$(c["observed_orders"]))")
        end
        if !c["conservation_pass"]
            push!(hard_fails, "conservation_pass failed for order case p=$p")
        end

        cc = run_advection_conservation(; p=p, n_elements=32)
        push!(cases, cc)
        if !cc["conservation_pass"]
            push!(hard_fails, "conservation residual too large for p=$p: $(cc["conservation_residual"])")
        end
    end

    overall = all(c -> c["pass"] === true, cases) && isempty(hard_fails)
    return cases, overall, hard_fails
end

# ---------------------------------------------------------------------------
# Milestone 2 — inviscid Burgers
# ---------------------------------------------------------------------------

"""Periodic square-wave IC for Burgers oscillation demo on [0,1]."""
function burgers_square_ic(x; u_left=1.0, u_right=0.0, x_disc=0.5)
    return x < x_disc ? u_left : u_right
end

"""
    run_burgers_conservation(; p, n_elements, t_final, cfl)

Periodic Burgers with discontinuous IC; check mass conservation.
"""
function run_burgers_conservation(;
    p::Int=3,
    n_elements::Int=32,
    t_final::Float64=0.15,
    cfl::Float64=0.2,
    method::AbstractCapturingMethod=NullCapturing(),
)
    t0 = time()
    eq = Burgers1D()
    ops = build_operators(p)
    mesh = Mesh1D(0.0, 1.0, n_elements; left_bc=PeriodicBC(), right_bc=PeriodicBC())
    state = allocate_state(mesh, ops, Val(1))
    set_initial_condition!(state, x -> burgers_square_ic(x))
    M0 = discrete_mass(state, 1)
    result = ssp_rk3!(state, eq, method, t_final; cfl=cfl)
    MT = discrete_mass(state, 1)
    cres = conservation_residual_relative(M0, MT)
    abs_res = conservation_residual_absolute(M0, MT)
    cpass = result.status == :ok && (cres <= 1e-10 || abs_res <= 1e-12)
    diverged = result.status != :ok
    nan_detected = diverged || has_nonfinite(state.u)

    return Dict{String,Any}(
        "name" => "burgers_conservation_p$(p)",
        "case_type" => "other",
        "equation" => "burgers",
        "p" => p,
        "capturing_method" => "null",
        "pass" => cpass && !diverged,
        "diverged" => diverged,
        "nan_detected" => nan_detected,
        "conservation_residual" => cres,
        "conservation_pass" => cpass,
        "conservation_metric" => "periodic_mass_change",
        "positivity_ok" => true,
        "wall_time_sec" => time() - t0,
        "metrics" => Dict{String,Any}(
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

Success for M2 means the run completes, conserves, and **exhibits** Gibbs-type
oscillations (overshoot above IC max or undershoot below IC min). That documents
the failure mode that shock-capturing methods later aim to fix.
"""
function run_burgers_oscillation(;
    p::Int=3,
    n_elements::Int=32,
    t_final::Float64=0.15,
    cfl::Float64=0.2,
    min_overshoot::Float64=0.02,  # require at least 2% of jump height
    method::AbstractCapturingMethod=NullCapturing(),
    u_left::Float64=1.0,
    u_right::Float64=0.0,
)
    t0 = time()
    eq = Burgers1D()
    ops = build_operators(p)
    mesh = Mesh1D(0.0, 1.0, n_elements; left_bc=PeriodicBC(), right_bc=PeriodicBC())
    state = allocate_state(mesh, ops, Val(1))
    set_initial_condition!(state, x -> burgers_square_ic(x; u_left=u_left, u_right=u_right))

    u0_min, u0_max = solution_extrema(state, 1)
    M0 = discrete_mass(state, 1)

    result = ssp_rk3!(state, eq, method, t_final; cfl=cfl)
    MT = discrete_mass(state, 1)
    u_min, u_max = solution_extrema(state, 1)
    η_over, η_under, η = overshoot_metric(u_min, u_max, u0_min, u0_max)

    cres = conservation_residual_relative(M0, MT)
    abs_res = conservation_residual_absolute(M0, MT)
    cpass = result.status == :ok && (cres <= 1e-10 || abs_res <= 1e-12)
    diverged = result.status != :ok
    nan_detected = diverged || has_nonfinite(state.u)

    # M2 success: oscillations present (this is the expected HO failure mode)
    oscillations_present = η >= min_overshoot
    # Case "pass" means demo succeeded: ran, conserved, and showed oscillations
    case_pass = cpass && !diverged && !nan_detected && oscillations_present

    return Dict{String,Any}(
        "name" => "burgers_oscillation_p$(p)",
        "case_type" => "discontinuous",
        "equation" => "burgers",
        "p" => p,
        "capturing_method" => "null",
        "pass" => case_pass,
        "diverged" => diverged,
        "nan_detected" => nan_detected,
        "conservation_residual" => cres,
        "conservation_pass" => cpass,
        "conservation_metric" => "periodic_mass_change",
        "positivity_ok" => true,  # scalar Burgers; no density constraint
        "wall_time_sec" => time() - t0,
        "n_elements" => n_elements,
        "t_final" => t_final,
        "excess_dissipation" => nothing,  # no NullCapturing reference needed for M2 demo
        "shock_thickness" => nothing,
        "shock_thickness_unit" => "sp_spacings",
        "overshoot" => η,
        "metrics" => Dict{String,Any}(
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
            "note" => "M2 documents HO oscillatory failure of NullCapturing; oscillations are expected",
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
        cc = run_burgers_conservation(; p=p, n_elements=32, t_final=0.15)
        push!(cases, cc)
        if !cc["conservation_pass"]
            push!(hard_fails, "burgers conservation failed for p=$p: res=$(cc["conservation_residual"])")
        end

        co = run_burgers_oscillation(; p=p, n_elements=32, t_final=0.15)
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
# Milestone 3 — 1D Euler + BCs + smooth order
# ---------------------------------------------------------------------------

"""
Exact density-wave solution (periodic):
  ρ = 1 + A sin(2π(x - u0 t)),  u = u0,  p = p0  (constant)
Conserved: (ρ, ρu, E) with E = p/(γ-1) + ½ ρ u².
"""
function euler_density_wave_conserved(eq::Euler1D{T}, x, t; A=T(0.2), u0=T(1), p0=T(1)) where {T}
    ρ = one(T) + T(A) * sin(2π * (x - u0 * t))
    return primitives_to_conserved(eq, ρ, u0, p0)
end

function _euler_order_study_dt(p::Int, Ne_fine::Int, λ_max::Float64; cfl::Float64=0.2)
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
    p::Int=3,
    n_elements_list::AbstractVector{Int}=[8, 16, 32],
    t_final::Float64=1.0,
    cfl::Float64=0.2,
    order_tol::Float64=DEFAULT_ORDER_TOLERANCE,
    γ::Float64=1.4,
    method::AbstractCapturingMethod=NullCapturing(),
)
    t0 = time()
    eq = Euler1D(γ)
    ops = build_operators(p)
    formal = p + 1
    Ne_fine = maximum(n_elements_list)
    λ_max = 2.5
    dt_fixed = _euler_order_study_dt(p, Ne_fine, λ_max; cfl=cfl)

    mesh_sizes = Float64[]
    l2_errors = Float64[]
    diverged = false
    nan_detected = false
    pos_ok = true

    uexact = x -> euler_density_wave_conserved(eq, x, t_final)

    for Ne in n_elements_list
        mesh = Mesh1D(0.0, 1.0, Ne; left_bc=PeriodicBC(), right_bc=PeriodicBC())
        state = allocate_state(mesh, ops, Val(3))
        set_initial_condition!(state, x -> euler_density_wave_conserved(eq, x, 0.0))
        result = ssp_rk3!(state, eq, method, t_final; dt=dt_fixed)
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
    opass = !diverged && order_pass(obs; formal_order=formal, tol=order_tol)

    Ne = n_elements_list[end]
    mesh = Mesh1D(0.0, 1.0, Ne; left_bc=PeriodicBC(), right_bc=PeriodicBC())
    state = allocate_state(mesh, ops, Val(3))
    set_initial_condition!(state, x -> euler_density_wave_conserved(eq, x, 0.0))
    M0 = [discrete_mass(state, c) for c in 1:3]
    result = ssp_rk3!(state, eq, method, t_final; cfl=cfl)
    MT = [discrete_mass(state, c) for c in 1:3]
    cres_comp = [conservation_residual_relative(M0[c], MT[c]) for c in 1:3]
    abs_comp = [conservation_residual_absolute(M0[c], MT[c]) for c in 1:3]
    cres = maximum(cres_comp)
    cpass =
        result.status == :ok &&
        all(cres_comp[c] <= 1e-10 || abs_comp[c] <= 1e-11 for c in 1:3)
    pos_ok = pos_ok && positivity_ok(eq, state)

    case_pass = opass && cpass && !diverged && !nan_detected && pos_ok

    return Dict{String,Any}(
        "name" => "euler_smooth_order_p$(p)",
        "case_type" => "smooth_order",
        "equation" => "euler1d",
        "p" => p,
        "capturing_method" => "null",
        "pass" => case_pass,
        "diverged" => diverged,
        "nan_detected" => nan_detected,
        "conservation_residual" => cres,
        "conservation_pass" => cpass,
        "conservation_metric" => "periodic_mass_change",
        "positivity_ok" => pos_ok,
        "wall_time_sec" => time() - t0,
        "metrics" => Dict{String,Any}(
            "γ" => γ,
            "t_final" => t_final,
            "cfl" => cfl,
            "dt_fixed" => dt_fixed,
            "conservation_residual_components" => cres_comp,
            "mass_abs_change_components" => abs_comp,
            "n_elements_list" => collect(n_elements_list),
        ),
        "mesh_sizes" => mesh_sizes,
        "l2_errors" => l2_errors,
        "observed_orders" => obs,
        "formal_order" => formal,
        "order_pass" => opass,
        "order_tolerance" => order_tol,
    )
end

"""
    run_euler_conservation(; p, n_elements, t_final, cfl)

Periodic Euler density wave: mass/momentum/energy conservation.
"""
function run_euler_conservation(;
    p::Int=3,
    n_elements::Int=32,
    t_final::Float64=1.0,
    cfl::Float64=0.2,
    γ::Float64=1.4,
    method::AbstractCapturingMethod=NullCapturing(),
)
    t0 = time()
    eq = Euler1D(γ)
    ops = build_operators(p)
    mesh = Mesh1D(0.0, 1.0, n_elements; left_bc=PeriodicBC(), right_bc=PeriodicBC())
    state = allocate_state(mesh, ops, Val(3))
    set_initial_condition!(state, x -> euler_density_wave_conserved(eq, x, 0.0))
    M0 = [discrete_mass(state, c) for c in 1:3]
    result = ssp_rk3!(state, eq, method, t_final; cfl=cfl)
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

    return Dict{String,Any}(
        "name" => "euler_conservation_p$(p)",
        "case_type" => "other",
        "equation" => "euler1d",
        "p" => p,
        "capturing_method" => "null",
        "pass" => cpass && !diverged && pos,
        "diverged" => diverged,
        "nan_detected" => nan_detected,
        "conservation_residual" => cres,
        "conservation_pass" => cpass,
        "conservation_metric" => "periodic_mass_change",
        "positivity_ok" => pos,
        "wall_time_sec" => time() - t0,
        "metrics" => Dict{String,Any}(
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
    p::Int=3,
    n_elements::Int=16,
    t_final::Float64=0.2,
    cfl::Float64=0.2,
    γ::Float64=1.4,
)
    t0 = time()
    eq = Euler1D(γ)
    ops = build_operators(p)
    mesh = Mesh1D(
        0.0,
        1.0,
        n_elements;
        left_bc=TransmissiveBC(),
        right_bc=TransmissiveBC(),
    )
    state = allocate_state(mesh, ops, Val(3))
    U0 = primitives_to_conserved(eq, 1.0, 0.5, 1.0)
    set_initial_condition!(state, _ -> U0)
    result = ssp_rk3!(state, eq, NullCapturing(), t_final; cfl=cfl)
    dev = 0.0
    for c in 1:3
        for e in 1:n_elements, j in 1:size(state.u, 1)
            dev = max(dev, abs(state.u[j, e, c] - U0[c]))
        end
    end
    ok = result.status == :ok && dev < 1e-10 && positivity_ok(eq, state)

    return Dict{String,Any}(
        "name" => "bc_transmissive_freestream_p$(p)",
        "case_type" => "other",
        "equation" => "euler1d",
        "p" => p,
        "capturing_method" => "null",
        "pass" => ok,
        "diverged" => result.status != :ok,
        "nan_detected" => has_nonfinite(state.u),
        "conservation_residual" => dev,
        "conservation_pass" => dev < 1e-10,
        "conservation_metric" => "none",
        "positivity_ok" => positivity_ok(eq, state),
        "wall_time_sec" => time() - t0,
        "metrics" => Dict{String,Any}(
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
    p::Int=3,
    n_elements::Int=16,
    t_final::Float64=0.2,
    cfl::Float64=0.2,
    γ::Float64=1.4,
)
    t0 = time()
    eq = Euler1D(γ)
    ops = build_operators(p)
    U0 = primitives_to_conserved(eq, 1.0, 0.3, 1.0)
    mesh = Mesh1D(
        0.0,
        1.0,
        n_elements;
        left_bc=DirichletBC(t -> U0),
        right_bc=DirichletBC(t -> U0),
    )
    state = allocate_state(mesh, ops, Val(3))
    set_initial_condition!(state, _ -> U0)
    result = ssp_rk3!(state, eq, NullCapturing(), t_final; cfl=cfl)
    dev = 0.0
    for c in 1:3
        for e in 1:n_elements, j in 1:size(state.u, 1)
            dev = max(dev, abs(state.u[j, e, c] - U0[c]))
        end
    end
    ok = result.status == :ok && dev < 1e-10 && positivity_ok(eq, state)

    return Dict{String,Any}(
        "name" => "bc_dirichlet_freestream_p$(p)",
        "case_type" => "other",
        "equation" => "euler1d",
        "p" => p,
        "capturing_method" => "null",
        "pass" => ok,
        "diverged" => result.status != :ok,
        "nan_detected" => has_nonfinite(state.u),
        "conservation_residual" => dev,
        "conservation_pass" => dev < 1e-10,
        "conservation_metric" => "none",
        "positivity_ok" => positivity_ok(eq, state),
        "wall_time_sec" => time() - t0,
        "metrics" => Dict{String,Any}(
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
        c = run_euler_smooth_order(; p=p, n_elements_list=nlist, cfl=0.1)
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

        cc = run_euler_conservation(; p=p, n_elements=32)
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
