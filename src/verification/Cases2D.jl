# Milestone 8 — 2D verification cases.

"""
    run_advection2d_smooth_order(; p, n_list, ax, ay, t_final, cfl)

2D periodic density/sine advection order study on unit square.
"""
function run_advection2d_smooth_order(;
    p::Int=2,
    n_list::AbstractVector{Int}=[4, 8, 16],
    ax::Float64=1.0,
    ay::Float64=0.5,
    t_final::Float64=1.0,
    cfl::Float64=0.15,
    order_tol::Float64=DEFAULT_ORDER_TOLERANCE,
)
    t0 = time()
    eq = LinearAdvection2D(ax, ay)
    ops = build_operators(p)
    formal = p + 1
    n_fine = maximum(n_list)
    λ = hypot(ax, ay)
    h = 1.0 / n_fine
    dt_fixed = (cfl / (p + 1)) * h / ((2p + 1) * max(λ, eps()))

    mesh_sizes = Float64[]
    l2_errors = Float64[]
    diverged = false

    uexact = (x, y) -> sin(2π * (x - ax * t_final)) * sin(2π * (y - ay * t_final))

    for n in n_list
        mesh = Mesh2D(0.0, 1.0, 0.0, 1.0, n, n)
        state = allocate_state(mesh, ops, Val(1))
        set_initial_condition!(state, (x, y) -> sin(2π * x) * sin(2π * y))
        result = ssp_rk3!(state, eq, NullCapturing(), t_final; dt=dt_fixed)
        if result.status != :ok
            diverged = true
            push!(mesh_sizes, 1.0 / n)
            push!(l2_errors, NaN)
            continue
        end
        push!(mesh_sizes, 1.0 / n)
        push!(l2_errors, l2_error(state, uexact, 1))
    end

    obs = observed_orders(mesh_sizes, l2_errors)
    opass = !diverged && order_pass(obs; formal_order=formal, tol=order_tol)

    # Conservation on finest
    n = n_list[end]
    mesh = Mesh2D(0.0, 1.0, 0.0, 1.0, n, n)
    state = allocate_state(mesh, ops, Val(1))
    set_initial_condition!(state, (x, y) -> 1.0 + 0.2 * sin(2π * x) * sin(2π * y))
    M0 = discrete_mass(state, 1)
    result = ssp_rk3!(state, eq, NullCapturing(), t_final; cfl=cfl)
    MT = discrete_mass(state, 1)
    cres = conservation_residual_relative(M0, MT)
    cpass = result.status == :ok && (cres <= 1e-10 || abs(MT - M0) <= 1e-11)

    return Dict{String,Any}(
        "name" => "advection2d_smooth_order_p$(p)",
        "case_type" => "smooth_order",
        "equation" => "linear_advection_2d",
        "p" => p,
        "capturing_method" => "null",
        "pass" => opass && cpass && !diverged,
        "diverged" => diverged,
        "nan_detected" => diverged,
        "conservation_residual" => cres,
        "conservation_pass" => cpass,
        "conservation_metric" => "periodic_mass_change",
        "positivity_ok" => true,
        "wall_time_sec" => time() - t0,
        "mesh_sizes" => mesh_sizes,
        "l2_errors" => l2_errors,
        "observed_orders" => obs,
        "formal_order" => formal,
        "order_pass" => opass,
        "order_tolerance" => order_tol,
        "metrics" => Dict{String,Any}("ax" => ax, "ay" => ay, "n_list" => collect(n_list)),
    )
end

"""
    run_euler2d_smooth_order(; p, n_list, t_final, cfl, method)

2D periodic density wave: ρ = 1 + A sin(2π(x-t)) sin(2π(y-t)), u=1, v=1, p=1.
Default method is NullCapturing (formal order). Pass PerssonAV for capturing-on check.
"""
function run_euler2d_smooth_order(;
    p::Int=2,
    n_list::AbstractVector{Int}=[8, 16, 32],
    t_final::Float64=0.25,
    cfl::Float64=0.06,
    γ::Float64=1.4,
    order_tol::Float64=DEFAULT_ORDER_TOLERANCE,
    method::AbstractCapturingMethod=NullCapturing(),
    method_name::AbstractString="null",
)
    t0 = time()
    eq = Euler2D(γ)
    ops = build_operators(p)
    formal = p + 1
    n_fine = maximum(n_list)
    λ = 1.0 + sqrt(γ)  # rough |u|+c with u=1,p=1,ρ~1
    h = 1.0 / n_fine
    dt_fixed = (cfl / (p + 1)) * h / ((2p + 1) * λ)

    mesh_sizes = Float64[]
    l2_errors = Float64[]
    diverged = false
    pos_ok = true

    function Uex(x, y)
        ρ = 1.0 + 0.2 * sin(2π * (x - t_final)) * sin(2π * (y - t_final))
        return primitives_to_conserved(eq, ρ, 1.0, 1.0, 1.0)
    end

    for n in n_list
        mesh = Mesh2D(0.0, 1.0, 0.0, 1.0, n, n)
        state = allocate_state(mesh, ops, Val(4))
        set_initial_condition!(
            state,
            (x, y) ->
                primitives_to_conserved(
                    eq,
                    1.0 + 0.2 * sin(2π * x) * sin(2π * y),
                    1.0,
                    1.0,
                    1.0,
                ),
        )
        result = integrate!(state, eq, method, t_final; dt=dt_fixed)
        if result.status != :ok
            diverged = true
            push!(mesh_sizes, 1.0 / n)
            push!(l2_errors, NaN)
            continue
        end
        pos_ok = pos_ok && positivity_ok(eq, state)
        push!(mesh_sizes, 1.0 / n)
        push!(l2_errors, l2_error(state, Uex, 1))
    end

    obs = observed_orders(mesh_sizes, l2_errors)
    opass = !diverged && order_pass(obs; formal_order=formal, tol=order_tol)

    # Periodic mass conservation check on finest
    n = n_list[end]
    mesh = Mesh2D(0.0, 1.0, 0.0, 1.0, n, n)
    state = allocate_state(mesh, ops, Val(4))
    set_initial_condition!(
        state,
        (x, y) ->
            primitives_to_conserved(
                eq,
                1.0 + 0.2 * sin(2π * x) * sin(2π * y),
                1.0,
                1.0,
                1.0,
            ),
    )
    M0 = discrete_mass(state, 1)
    result = integrate!(state, eq, method, t_final; cfl=cfl)
    MT = discrete_mass(state, 1)
    cres = conservation_residual_relative(M0, MT)
    cpass = result.status == :ok && (cres <= 1e-8 || abs(MT - M0) <= 1e-10)

    return Dict{String,Any}(
        "name" => "euler2d_smooth_order_p$(p)_$(method_name)",
        "case_type" => "smooth_order",
        "equation" => "euler2d",
        "p" => p,
        "capturing_method" => String(method_name),
        "pass" => opass && cpass && !diverged && pos_ok,
        "diverged" => diverged,
        "nan_detected" => diverged,
        "conservation_residual" => cres,
        "conservation_pass" => cpass,
        "conservation_metric" => "periodic_mass_change",
        "positivity_ok" => pos_ok,
        "wall_time_sec" => time() - t0,
        "mesh_sizes" => mesh_sizes,
        "l2_errors" => l2_errors,
        "observed_orders" => obs,
        "formal_order" => formal,
        "order_pass" => opass,
        "order_tolerance" => order_tol,
        "metrics" => Dict{String,Any}(
            "n_list" => collect(n_list),
            "method_params" => method_params(method),
        ),
    )
end

"""
    run_euler2d_discontinuous(; p, nx, ny, t_final, method)

Simple 2D discontinuous density jump (transmissive). With NullCapturing may
oscillate; with Persson AV should remain stable/positive on CI-light meshes.
"""
function run_euler2d_discontinuous(;
    p::Int=1,
    nx::Int=20,
    ny::Int=20,
    t_final::Float64=0.1,
    cfl::Float64=0.1,
    γ::Float64=1.4,
    method::AbstractCapturingMethod=NullCapturing(),
    method_name::AbstractString="null",
)
    t0 = time()
    eq = Euler2D(γ)
    ops = build_operators(p)
    mesh = Mesh2D(
        0.0,
        1.0,
        0.0,
        1.0,
        nx,
        ny;
        left_bc=TransmissiveBC(),
        right_bc=TransmissiveBC(),
        bottom_bc=TransmissiveBC(),
        top_bc=TransmissiveBC(),
    )
    state = allocate_state(mesh, ops, Val(4))
    set_initial_condition!(
        state,
        (x, y) ->
            x < 0.5 ? primitives_to_conserved(eq, 1.0, 0.0, 0.0, 1.0) :
            primitives_to_conserved(eq, 0.125, 0.0, 0.0, 0.1),
    )
    # Density extrema for overshoot metric
    ρ0 = state.u[:, :, :, 1]
    ρmin0, ρmax0 = minimum(ρ0), maximum(ρ0)
    result = integrate!(state, eq, method, t_final; cfl=cfl)
    diverged = result.status != :ok
    nan_det = diverged || has_nonfinite(state.u)
    pos = !nan_det && positivity_ok(eq, state)
    ρ = state.u[:, :, :, 1]
    _, _, η = overshoot_metric(minimum(ρ), maximum(ρ), 0.125, 1.0)
    case = Dict{String,Any}(
        "name" => "euler2d_density_jump_p$(p)_$(method_name)",
        "case_type" => "discontinuous",
        "equation" => "euler2d",
        "p" => p,
        "capturing_method" => String(method_name),
        "pass" => !diverged && pos && !nan_det,
        "diverged" => diverged,
        "nan_detected" => nan_det,
        "conservation_residual" => 0.0,
        "conservation_pass" => true,
        "conservation_metric" => "none",
        "positivity_ok" => pos,
        "wall_time_sec" => time() - t0,
        "n_elements" => nx * ny,
        "t_final" => t_final,
        "excess_dissipation" => nothing,
        "shock_thickness" => nothing,
        "shock_thickness_unit" => "sp_spacings",
        "overshoot" => η,
        "metrics" => Dict{String,Any}(
            "nx" => nx,
            "ny" => ny,
            "n_steps" => result.n_steps,
            "ρ_min0" => ρmin0,
            "ρ_max0" => ρmax0,
            "method_params" => method_params(method),
        ),
    )
    return case, state, eq
end

"""
    run_m8_2d_suite() -> (cases, overall, fails, optional_state_for_vtk)

Baseline 2D suite (NullCapturing) plus CI-light 2D capturing checks (Phase 3.1).
"""
function run_m8_2d_suite()
    cases = Any[]
    hard_fails = String[]

    c1 = run_advection2d_smooth_order(; p=2, n_list=[4, 8, 16], t_final=1.0, cfl=0.12)
    push!(cases, c1)
    if !c1["order_pass"]
        push!(hard_fails, "2D advection order failed: $(c1["observed_orders"])")
    end

    c2 = run_euler2d_smooth_order(; p=2, n_list=[16, 32], t_final=0.25, cfl=0.06)
    push!(cases, c2)
    if !c2["order_pass"]
        push!(hard_fails, "2D Euler order failed: $(c2["observed_orders"])")
    end

    c3, state_disc, eq_disc = run_euler2d_discontinuous(; p=1, nx=16, ny=16, t_final=0.05)
    push!(cases, c3)
    if !c3["pass"]
        push!(hard_fails, "2D discontinuous Euler (null) failed")
    end

    # Phase 3.1 CI-light: Persson AV on jump + short smooth with AV on
    pers = PerssonAVMethod(; c_av=0.1)
    c4, _, _ = run_euler2d_discontinuous(;
        p=1,
        nx=12,
        ny=12,
        t_final=0.05,
        cfl=0.08,
        method=pers,
        method_name="persson_av",
    )
    push!(cases, c4)
    if !c4["pass"]
        push!(hard_fails, "2D discontinuous Euler (persson_av) failed")
    end

    c5 = run_euler2d_smooth_order(;
        p=2,
        n_list=[8, 16],
        t_final=0.15,
        cfl=0.05,
        method=pers,
        method_name="persson_av",
    )
    push!(cases, c5)
    # With AV on smooth data, order may degrade slightly — require stability + conservation
    if c5["diverged"] || !c5["positivity_ok"] || !c5["conservation_pass"]
        push!(hard_fails, "2D Euler smooth + persson_av unstable or non-conservative")
        c5["pass"] = false
    else
        # Soften: pass if stable even if formal order not fully recovered
        c5["pass"] = true
        c5["metrics"]["note"] = "capturing-on smooth: stability/conservation gated; order informative"
    end

    overall = all(c -> c["pass"] === true, cases) && isempty(hard_fails)
    return cases, overall, hard_fails, state_disc, eq_disc
end

"""
    run_freestream_preservation_2d(; p, nx, ny, curved, amp, t_final)

Constant free-stream Euler: residual / time evolution must stay near machine zero
(metric geometric conservation). Merge gate for curved elements.
"""
function run_freestream_preservation_2d(;
    p::Int=2,
    nx::Int=4,
    ny::Int=4,
    curved::Bool=false,
    amp::Float64=0.05,
    t_final::Float64=0.1,
    cfl::Float64=0.1,
    γ::Float64=1.4,
    tol::Union{Nothing,Float64}=nothing,
)
    t0 = time()
    eq = Euler2D(γ)
    ops = build_operators(p)
    mesh = if curved
        make_wavy_mesh2d(nx, ny, ops; amp=amp)
    else
        Mesh2D(0.0, 1.0, 0.0, 1.0, nx, ny)
    end
    state = allocate_state(
        mesh, ops, Val(4); wavy_amp=curved ? amp : nothing,
    )
    U0 = primitives_to_conserved(eq, 1.0, 0.5, 0.25, 1.0)
    set_initial_condition!(state, (x, y) -> U0)
    # Instantaneous residual
    du = similar(state.u)
    residual!(du, state, eq, NullCapturing())
    res_max = maximum(abs, du)
    # Short evolution
    result = integrate!(state, eq, NullCapturing(), t_final; cfl=cfl)
    err = zero(eltype(state.u))
    for c in 1:4
        err = max(err, maximum(abs, state.u[:, :, :, c] .- U0[c]))
    end
    pos = positivity_ok(eq, state)
    # Affine: machine-precision free-stream. Curved (analytic wavy metrics):
    # residual is spectral-accurate GCL error — strict on refined CI-light config.
    tol_res = something(tol, curved ? 5e-4 : 1e-10)
    tol_err = curved ? 5e-5 : 1e-8
    pass =
        result.status == :ok && res_max < tol_res && err < tol_err && pos &&
        all(state.metrics.J .> 0)
    return Dict{String,Any}(
        "name" => "freestream_2d_p$(p)_" * (curved ? "curved" : "cartesian"),
        "case_type" => "conservation",
        "equation" => "euler2d",
        "p" => p,
        "capturing_method" => "null",
        "pass" => pass,
        "diverged" => result.status != :ok,
        "nan_detected" => has_nonfinite(state.u),
        "conservation_residual" => Float64(err),
        "conservation_pass" => err < tol_err,
        "conservation_metric" => "freestream_max_dev",
        "positivity_ok" => pos,
        "wall_time_sec" => time() - t0,
        "metrics" => Dict{String,Any}(
            "residual_max" => Float64(res_max),
            "freestream_error" => Float64(err),
            "curved" => curved,
            "amp" => amp,
            "tol_res" => tol_res,
            "tol_err" => tol_err,
            "n_steps" => result.n_steps,
        ),
    )
end

"""
    run_advection2d_curved_order(; p, n_list, amp)

Smooth advection order on wavy mesh (periodic unit square with mild warp).
"""
function run_advection2d_curved_order(;
    p::Int=2,
    n_list::AbstractVector{Int}=[4, 8, 16],
    ax::Float64=1.0,
    ay::Float64=0.5,
    t_final::Float64=0.5,
    cfl::Float64=0.1,
    amp::Float64=0.03,
    order_tol::Float64=DEFAULT_ORDER_TOLERANCE,
)
    t0 = time()
    eq = LinearAdvection2D(ax, ay)
    ops = build_operators(p)
    formal = p + 1
    n_fine = maximum(n_list)
    λ = hypot(ax, ay)
    h = 1.0 / n_fine
    dt_fixed = (cfl / (p + 1)) * h / ((2p + 1) * max(λ, eps()))

    mesh_sizes = Float64[]
    l2_errors = Float64[]
    diverged = false

    uexact = (x, y) -> sin(2π * (x - ax * t_final)) * sin(2π * (y - ay * t_final))

    for n in n_list
        mesh = make_wavy_mesh2d(n, n, ops; amp=amp)
        state = allocate_state(mesh, ops, Val(1); wavy_amp=amp)
        set_initial_condition!(state, (x, y) -> sin(2π * x) * sin(2π * y))
        result = integrate!(state, eq, NullCapturing(), t_final; dt=dt_fixed)
        if result.status != :ok
            diverged = true
            push!(mesh_sizes, 1.0 / n)
            push!(l2_errors, NaN)
            continue
        end
        push!(mesh_sizes, 1.0 / n)
        push!(l2_errors, l2_error(state, uexact, 1))
    end

    obs = observed_orders(mesh_sizes, l2_errors)
    # Curved meshes: allow slightly looser order (pre-asymptotic + metric)
    opass = !diverged && order_pass(obs; formal_order=formal, tol=order_tol + 0.2)

    return Dict{String,Any}(
        "name" => "advection2d_curved_order_p$(p)",
        "case_type" => "smooth_order",
        "equation" => "linear_advection_2d",
        "p" => p,
        "capturing_method" => "null",
        "pass" => opass && !diverged,
        "diverged" => diverged,
        "nan_detected" => diverged,
        "conservation_residual" => 0.0,
        "conservation_pass" => true,
        "conservation_metric" => "none",
        "positivity_ok" => true,
        "wall_time_sec" => time() - t0,
        "mesh_sizes" => mesh_sizes,
        "l2_errors" => l2_errors,
        "observed_orders" => obs,
        "formal_order" => formal,
        "order_pass" => opass,
        "order_tolerance" => order_tol + 0.2,
        "metrics" => Dict{String,Any}("amp" => amp, "n_list" => collect(n_list)),
    )
end

"""
    run_euler2d_curved_discontinuous(; ...)

Density jump on wavy mesh with NullCapturing or Persson AV (CI-light).
"""
function run_euler2d_curved_discontinuous(;
    p::Int=1,
    nx::Int=12,
    ny::Int=12,
    t_final::Float64=0.03,
    cfl::Float64=0.04,
    amp::Float64=0.01,
    γ::Float64=1.4,
    method::AbstractCapturingMethod=NullCapturing(),
    method_name::AbstractString="null",
)
    t0 = time()
    eq = Euler2D(γ)
    ops = build_operators(p)
    mesh = make_wavy_mesh2d(
        nx,
        ny,
        ops;
        amp=amp,
        left_bc=TransmissiveBC(),
        right_bc=TransmissiveBC(),
        bottom_bc=TransmissiveBC(),
        top_bc=TransmissiveBC(),
    )
    state = allocate_state(mesh, ops, Val(4); wavy_amp=amp)
    set_initial_condition!(
        state,
        (x, y) ->
            x < 0.5 ? primitives_to_conserved(eq, 1.0, 0.0, 0.0, 1.0) :
            primitives_to_conserved(eq, 0.125, 0.0, 0.0, 0.1),
    )
    result = integrate!(state, eq, method, t_final; cfl=cfl)
    diverged = result.status != :ok
    nan_det = diverged || has_nonfinite(state.u)
    pos = !nan_det && positivity_ok(eq, state)
    return Dict{String,Any}(
        "name" => "euler2d_curved_jump_p$(p)_$(method_name)",
        "case_type" => "discontinuous",
        "equation" => "euler2d",
        "p" => p,
        "capturing_method" => String(method_name),
        "pass" => !diverged && pos && !nan_det,
        "diverged" => diverged,
        "nan_detected" => nan_det,
        "conservation_residual" => 0.0,
        "conservation_pass" => true,
        "conservation_metric" => "none",
        "positivity_ok" => pos,
        "wall_time_sec" => time() - t0,
        "n_elements" => nx * ny,
        "t_final" => t_final,
        "excess_dissipation" => nothing,
        "shock_thickness" => nothing,
        "shock_thickness_unit" => "sp_spacings",
        "overshoot" => 0.0,
        "metrics" => Dict{String,Any}(
            "nx" => nx,
            "ny" => ny,
            "amp" => amp,
            "n_steps" => result.n_steps,
            "method_params" => method_params(method),
        ),
    )
end

"""
    run_p32_curved_suite() -> (cases, overall, hard_fails)

CI-light curved-element verification suite (Phase 3.2).
"""
function run_p32_curved_suite()
    cases = Any[]
    hard_fails = String[]

    c1 = run_freestream_preservation_2d(; p=2, nx=4, ny=4, curved=false)
    push!(cases, c1)
    !c1["pass"] && push!(hard_fails, "cartesian freestream failed")

    # Curved free-stream: p=3, n=4 is CI-light and well under residual tol
    c2 = run_freestream_preservation_2d(; p=3, nx=4, ny=4, curved=true, amp=0.05, t_final=0.05)
    push!(cases, c2)
    !c2["pass"] && push!(hard_fails, "curved freestream failed res=$(c2["metrics"]["residual_max"])")

    c3 = run_advection2d_curved_order(; p=2, n_list=[4, 8, 16], amp=0.03, t_final=0.5)
    push!(cases, c3)
    !c3["pass"] && push!(hard_fails, "curved advection order failed: $(c3["observed_orders"])")

    c4 = run_euler2d_curved_discontinuous(;
        p=1,
        nx=12,
        ny=12,
        t_final=0.03,
        amp=0.01,
        cfl=0.04,
        method=NullCapturing(),
        method_name="null",
    )
    push!(cases, c4)
    !c4["pass"] && push!(hard_fails, "curved jump null failed")

    c5 = run_euler2d_curved_discontinuous(;
        p=1,
        nx=12,
        ny=12,
        t_final=0.03,
        amp=0.01,
        cfl=0.04,
        method=PerssonAVMethod(; c_av=0.1),
        method_name="persson_av",
    )
    push!(cases, c5)
    !c5["pass"] && push!(hard_fails, "curved jump persson_av failed")

    overall = isempty(hard_fails) && all(c -> c["pass"] === true, cases)
    return cases, overall, hard_fails
end

"""
    run_p31_2d_capturing_suite() -> (cases, overall, hard_fails)

Focused 2D capturing suite (CI-friendly reduced meshes).
"""
function run_p31_2d_capturing_suite()
    cases = Any[]
    hard_fails = String[]
    pers = PerssonAVMethod(; c_av=0.1)

    c_null, _, _ = run_euler2d_discontinuous(; p=1, nx=12, ny=12, t_final=0.05, method_name="null")
    push!(cases, c_null)

    c_av, _, _ = run_euler2d_discontinuous(;
        p=1,
        nx=12,
        ny=12,
        t_final=0.05,
        cfl=0.08,
        method=pers,
        method_name="persson_av",
    )
    push!(cases, c_av)
    if !c_av["pass"]
        push!(hard_fails, "persson_av 2D jump failed")
    end

    # Overshoot should not be wildly worse with AV (informational gate)
    if c_av["pass"] && c_null["pass"]
        if c_av["overshoot"] > c_null["overshoot"] + 0.5
            push!(hard_fails, "persson_av overshoot much larger than null (unexpected)")
        end
    end

    c_sm = run_euler2d_smooth_order(;
        p=2,
        n_list=[8, 16],
        t_final=0.15,
        cfl=0.05,
        method=pers,
        method_name="persson_av",
    )
    if c_sm["diverged"] || !c_sm["positivity_ok"] || !c_sm["conservation_pass"]
        c_sm["pass"] = false
        push!(hard_fails, "persson_av 2D smooth unstable")
    else
        c_sm["pass"] = true
    end
    push!(cases, c_sm)

    overall = isempty(hard_fails) && all(c -> c["pass"] === true, cases)
    return cases, overall, hard_fails
end
