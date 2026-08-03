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
