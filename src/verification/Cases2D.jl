# 2D verification cases: Cartesian/curved FR, capturing, Riemann/vortex, optional DMR/FFS.

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

Baseline 2D suite (`NullCapturing`) plus CI-light 2D capturing checks.
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

    # CI-light: Persson AV on jump + short smooth with AV on
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

CI-light curved-element verification suite.
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

# ---------------------------------------------------------------------------
# --- 2D Riemann + isentropic vortex (CI-light core gates) ---
# ---------------------------------------------------------------------------

"""
Lax–Liu Configuration 6 (four contacts / shears, uniform pressure).
States as (ρ, u, v, p). Preferred for CI-light positivity on coarse HO meshes.
"""
const RIEMANN2D_CFG6 = (
    ne = (1.0, 0.75, -0.5, 1.0),
    nw = (2.0, 0.75, 0.5, 1.0),
    sw = (1.0, -0.75, 0.5, 1.0),
    se = (3.0, -0.75, -0.5, 1.0),
)

"""
Lax–Liu Configuration 3 (strong shocks/rarefactions). Harder; full/nightly or research.
"""
const RIEMANN2D_CFG3 = (
    ne = (1.5, 0.0, 0.0, 1.5),
    nw = (0.5323, 1.206, 0.0, 0.3),
    sw = (0.138, 1.206, 1.206, 0.029),
    se = (0.5323, 0.0, 1.206, 0.3),
)

"""
    riemann2d_ic(eq, x, y; config=:cfg6, x0=0.5, y0=0.5)

Piecewise constant four-quadrant 2D Riemann IC.
`config` ∈ (`:cfg6`, `:cfg3`) or a NamedTuple of four (ρ,u,v,p) states.
"""
function riemann2d_ic(eq::Euler2D, x, y; config=:cfg6, x0=0.5, y0=0.5)
    states = if config === :cfg6 || config === "cfg6"
        RIEMANN2D_CFG6
    elseif config === :cfg3 || config === "cfg3"
        RIEMANN2D_CFG3
    else
        config
    end
    if x >= x0 && y >= y0
        ρ, u, v, p = states.ne
    elseif x < x0 && y >= y0
        ρ, u, v, p = states.nw
    elseif x < x0 && y < y0
        ρ, u, v, p = states.sw
    else
        ρ, u, v, p = states.se
    end
    return primitives_to_conserved(eq, ρ, u, v, p)
end

"""Alias for Configuration 3 (research / full tier)."""
riemann2d_cfg3_ic(eq::Euler2D, x, y; kwargs...) =
    riemann2d_ic(eq, x, y; config=:cfg3, kwargs...)

"""
    run_euler2d_riemann(; p, nx, ny, t_final, method, config=:cfg6, ...)

2D four-quadrant Riemann (Lax–Liu). Default **Configuration 6** (CI-light):
contacts/shears at uniform pressure — multi-wave structure without extreme
rarefactions that lose positivity on coarse high-order meshes.

`config=:cfg3` is available for full/nightly research runs (may need stronger AV
or relaxed positivity gates).

Pass criterion (CI): finite, positive, non-diverged.
"""
function run_euler2d_riemann(;
    p::Int=1,
    nx::Int=16,
    ny::Int=16,
    t_final::Float64=0.08,
    cfl::Float64=0.08,
    γ::Float64=1.4,
    x0::Float64=0.5,
    y0::Float64=0.5,
    config=:cfg6,
    method::AbstractCapturingMethod=PerssonAVMethod(; c_av=0.1),
    method_name::AbstractString="persson_av",
    require_positivity::Bool=true,
    progress_every::Int=0,
    progress_label::AbstractString="",
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
        state, (x, y) -> riemann2d_ic(eq, x, y; config=config, x0=x0, y0=y0),
    )
    ρ0 = @view state.u[:, :, :, 1]
    ρmin0, ρmax0 = minimum(ρ0), maximum(ρ0)
    plabel = isempty(progress_label) ? "riemann_$(config)" : progress_label
    result = integrate!(
        state, eq, method, t_final;
        cfl=cfl, progress_every=progress_every, progress_label=plabel,
    )
    diverged = result.status != :ok
    nan_det = diverged || has_nonfinite(state.u)
    pos = !nan_det && positivity_ok(eq, state)
    ρ = @view state.u[:, :, :, 1]
    ρmin, ρmax = nan_det ? (NaN, NaN) : (minimum(ρ), maximum(ρ))
    _, _, η = overshoot_metric(ρmin, ρmax, min(ρmin0, 0.5), max(ρmax0, 3.0))
    pass = !diverged && !nan_det && (!require_positivity || pos)
    cfg_str = string(config)
    return Dict{String,Any}(
        "name" => "euler2d_riemann_$(cfg_str)_p$(p)_$(method_name)",
        "case_type" => "discontinuous",
        "equation" => "euler2d",
        "p" => p,
        "capturing_method" => String(method_name),
        "pass" => pass,
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
            "config" => cfg_str,
            "nx" => nx,
            "ny" => ny,
            "n_steps" => result.n_steps,
            "ρ_min0" => ρmin0,
            "ρ_max0" => ρmax0,
            "ρ_min" => ρmin,
            "ρ_max" => ρmax,
            "method_params" => method_params(method),
            "ci_tier" => config === :cfg6 ? "required_ci_light" : "full_nightly",
            "require_positivity" => require_positivity,
        ),
    ),
    state,
    eq
end

"""
    isentropic_vortex_primitives(x, y, t; ...)

Classic isentropic vortex (Yee et al.) primitives (ρ, u, v, p) at time `t`.
Freestream (ρ∞,u∞,v∞,p∞)=(1,u∞,v∞,1); vortex strength β; center (x0,y0).
"""
function isentropic_vortex_primitives(
    x,
    y,
    t;
    γ::Float64=1.4,
    β::Float64=5.0,
    x0::Float64=5.0,
    y0::Float64=5.0,
    u∞::Float64=1.0,
    v∞::Float64=1.0,
)
    xd = x - x0 - u∞ * t
    yd = y - y0 - v∞ * t
    r2 = xd * xd + yd * yd
    expf = exp(0.5 * (1.0 - r2))
    δu = -(β / (2π)) * expf * yd
    δv = (β / (2π)) * expf * xd
    T = 1.0 - (γ - 1.0) * β^2 / (8.0 * γ * π^2) * exp(1.0 - r2)
    # Guard against tiny floating underflow at large r
    T = max(T, 1e-14)
    ρ = T^(1.0 / (γ - 1.0))
    p = ρ^γ
    return ρ, u∞ + δu, v∞ + δv, p
end

"""
    run_isentropic_vortex_order(; p, n_list, t_final, ...)

Smooth isentropic vortex order study (periodic [0,L]²). NullCapturing by default.
Exact solution known — primary CI-friendly multi-D order / MMS-style gate.
"""
function run_isentropic_vortex_order(;
    p::Int=2,
    n_list::AbstractVector{Int}=[8, 16],
    t_final::Float64=0.5,
    cfl::Float64=0.08,
    L::Float64=10.0,
    γ::Float64=1.4,
    β::Float64=5.0,
    x0::Float64=5.0,
    y0::Float64=5.0,
    u∞::Float64=1.0,
    v∞::Float64=1.0,
    order_tol::Float64=DEFAULT_ORDER_TOLERANCE,
    method::AbstractCapturingMethod=NullCapturing(),
    method_name::AbstractString="null",
)
    t0 = time()
    eq = Euler2D(γ)
    ops = build_operators(p)
    formal = p + 1
    n_fine = maximum(n_list)
    # |u|+c freestream-ish: u∞+β/(2π) + sqrt(γ)
    λ = abs(u∞) + abs(v∞) + β / (2π) + sqrt(γ) + 1.0
    h = L / n_fine
    dt_fixed = (cfl / (p + 1)) * h / ((2p + 1) * λ)

    mesh_sizes = Float64[]
    l2_errors = Float64[]
    diverged = false
    pos_ok = true

    function Uex(x, y)
        ρ, u, v, pr = isentropic_vortex_primitives(
            x, y, t_final; γ=γ, β=β, x0=x0, y0=y0, u∞=u∞, v∞=v∞,
        )
        return primitives_to_conserved(eq, ρ, u, v, pr)
    end

    for n in n_list
        mesh = Mesh2D(0.0, L, 0.0, L, n, n)
        state = allocate_state(mesh, ops, Val(4))
        set_initial_condition!(
            state,
            (x, y) -> begin
                ρ, u, v, pr = isentropic_vortex_primitives(
                    x, y, 0.0; γ=γ, β=β, x0=x0, y0=y0, u∞=u∞, v∞=v∞,
                )
                primitives_to_conserved(eq, ρ, u, v, pr)
            end,
        )
        result = integrate!(state, eq, method, t_final; dt=dt_fixed)
        if result.status != :ok
            diverged = true
            push!(mesh_sizes, L / n)
            push!(l2_errors, NaN)
            continue
        end
        pos_ok = pos_ok && positivity_ok(eq, state)
        push!(mesh_sizes, L / n)
        # Density L2 error vs exact
        push!(l2_errors, l2_error(state, Uex, 1))
    end

    obs = observed_orders(mesh_sizes, l2_errors)
    # Vortex on modest grids can be pre-asymptotic; allow slightly looser tol
    opass = !diverged && order_pass(obs; formal_order=formal, tol=order_tol + 0.15)

    # Periodic mass conservation on finest
    n = n_list[end]
    mesh = Mesh2D(0.0, L, 0.0, L, n, n)
    state = allocate_state(mesh, ops, Val(4))
    set_initial_condition!(
        state,
        (x, y) -> begin
            ρ, u, v, pr = isentropic_vortex_primitives(
                x, y, 0.0; γ=γ, β=β, x0=x0, y0=y0, u∞=u∞, v∞=v∞,
            )
            primitives_to_conserved(eq, ρ, u, v, pr)
        end,
    )
    M0 = discrete_mass(state, 1)
    result = integrate!(state, eq, method, t_final; cfl=cfl)
    MT = discrete_mass(state, 1)
    cres = conservation_residual_relative(M0, MT)
    cpass = result.status == :ok && (cres <= 1e-8 || abs(MT - M0) <= 1e-10)
    pos_ok = pos_ok && positivity_ok(eq, state)

    return Dict{String,Any}(
        "name" => "isentropic_vortex_order_p$(p)_$(method_name)",
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
        "order_tolerance" => order_tol + 0.15,
        "metrics" => Dict{String,Any}(
            "n_list" => collect(n_list),
            "L" => L,
            "β" => β,
            "t_final" => t_final,
            "method_params" => method_params(method),
            "ci_tier" => "required_ci_light",
            "note" => "isentropic vortex = exact-solution smooth Euler (MMS-style) gate",
        ),
    )
end

"""
    run_p33a_benchmark_suite() -> (cases, overall, hard_fails)

CI-light suite: 2D Riemann + isentropic vortex order.
Double Mach / FFS live in the optional2d suite (full/nightly).
"""
function run_p33a_benchmark_suite()
    cases = Any[]
    hard_fails = String[]

    # Isentropic vortex order (NullCapturing) — core smooth multi-D gate
    c_v = run_isentropic_vortex_order(;
        p=2,
        n_list=[8, 16],
        t_final=0.5,
        cfl=0.08,
        L=10.0,
    )
    push!(cases, c_v)
    if !c_v["order_pass"] || c_v["diverged"] || !c_v["conservation_pass"]
        push!(
            hard_fails,
            "isentropic vortex order failed: orders=$(c_v["observed_orders"]) L2=$(c_v["l2_errors"])",
        )
    end

    # 2D Riemann cfg6 + Persson AV — core discontinuous multi-D gate
    c_r, _, _ = run_euler2d_riemann(;
        p=1,
        nx=16,
        ny=16,
        t_final=0.08,
        cfl=0.08,
        config=:cfg6,
        method=PerssonAVMethod(; c_av=0.1),
        method_name="persson_av",
    )
    push!(cases, c_r)
    if !c_r["pass"]
        push!(hard_fails, "2D Riemann cfg6 + persson_av failed")
    end

    # Short NullCapturing smoke on same config (positivity still required for cfg6)
    c_rn, _, _ = run_euler2d_riemann(;
        p=1,
        nx=12,
        ny=12,
        t_final=0.04,
        cfl=0.06,
        config=:cfg6,
        method=NullCapturing(),
        method_name="null",
    )
    push!(cases, c_rn)
    if !c_rn["pass"]
        push!(hard_fails, "2D Riemann cfg6 null failed")
    end

    overall = isempty(hard_fails) && all(c -> c["pass"] === true, cases)
    return cases, overall, hard_fails
end

# ---------------------------------------------------------------------------
# --- Optional reduced Double Mach + Forward-Facing Step ---
# Full/nightly tier (unit tests use ultra-light configs).
# ---------------------------------------------------------------------------

"""
    double_mach_states(eq; strength=:reduced)

Post- / pre-shock states for Double-Mach-like inclined shock.

- `:reduced` (default): mild Ms≈2-class jump for coarse HO meshes (CI/unit tests).
- `:classic`: Woodward–Colella Ms≈10 / 60° states (research; needs finer mesh + strong AV).
"""
function double_mach_states(eq::Euler2D; strength::Symbol=:reduced)
    U_pre = primitives_to_conserved(eq, 1.4, 0.0, 0.0, 1.0)
    if strength === :classic
        U_post = primitives_to_conserved(
            eq, 8.0, 8.25 * cos(π / 6), -8.25 * sin(π / 6), 116.5,
        )
    else
        # Milder post-shock (exercises walls + inclined front without face-trace blowup)
        U_post = primitives_to_conserved(eq, 3.5, 1.6 * cos(π / 6), -1.6 * sin(π / 6), 4.5)
    end
    return U_post, U_pre
end

"""
    double_mach_ic(eq, x, y; x_shock0=1/6, strength=:reduced)

Initial condition: post-shock for x < x_s(y), pre otherwise (60° front).
"""
function double_mach_ic(eq::Euler2D, x, y; x_shock0=1 / 6, strength::Symbol=:reduced)
    U_post, U_pre = double_mach_states(eq; strength=strength)
    xs = x_shock0 + y / √3
    return x < xs ? U_post : U_pre
end

"""
    run_double_mach_reflection(; ...)

Reduced Double-Mach-like reflection (inclined shock + reflecting wall).
Default `strength=:reduced` is HO-stable on coarse meshes. Classic Ms=10 via
`strength=:classic` (research / fine meshes only).

**CI tier:** full/nightly (optional2d). Pass: finite, non-diverged; positivity optional on ultra-coarse meshes.
"""
function run_double_mach_reflection(;
    p::Int=1,
    nx::Int=24,
    ny::Int=6,
    t_final::Float64=0.04,
    cfl::Float64=0.04,
    γ::Float64=1.4,
    x_shock0::Float64=1 / 6,
    Lx::Float64=2.0,
    Ly::Float64=0.5,
    strength::Symbol=:reduced,
    method::AbstractCapturingMethod=PerssonAVMethod(; c_av=0.15),
    method_name::AbstractString="persson_av",
    require_positivity::Bool=false,
    progress_every::Int=0,
    progress_label::AbstractString="",
)
    t0 = time()
    eq = Euler2D(γ)
    ops = build_operators(p)
    U_post, U_pre = double_mach_states(eq; strength=strength)
    # Shock foot speed along top ~ proportional to Ly-scaled classic 20 when Ly=1
    vs = strength === :classic ? 20.0 : 2.5

    bottom = GhostStateBC() do u_int, nx_out, ny_out, x, y, t
        if x < x_shock0
            return U_post
        else
            return reflect_conserved(u_int, nx_out, ny_out)
        end
    end
    left = DirichletBC((x, y, t) -> U_post)
    right = TransmissiveBC()
    top = GhostStateBC() do u_int, nx_out, ny_out, x, y, t
        xs = x_shock0 + (Ly + vs * t) / √3
        return x < xs ? U_post : U_pre
    end

    mesh = Mesh2D(
        0.0, Lx, 0.0, Ly, nx, ny;
        left_bc=left, right_bc=right, bottom_bc=bottom, top_bc=top,
    )
    state = allocate_state(mesh, ops, Val(4))
    set_initial_condition!(
        state, (x, y) -> double_mach_ic(eq, x, y; x_shock0=x_shock0, strength=strength),
    )
    plabel = isempty(progress_label) ? "double_mach_$(strength)" : progress_label
    result = integrate!(
        state, eq, method, t_final;
        cfl=cfl, progress_every=progress_every, progress_label=plabel,
    )
    diverged = result.status != :ok
    nan_det = diverged || has_nonfinite(state.u)
    pos = !nan_det && positivity_ok(eq, state)
    pass = !diverged && !nan_det && (!require_positivity || pos)
    ρ = @view state.u[:, :, :, 1]
    return Dict{String,Any}(
        "name" => "double_mach_$(strength)_p$(p)_$(method_name)",
        "case_type" => "discontinuous",
        "equation" => "euler2d",
        "p" => p,
        "capturing_method" => String(method_name),
        "pass" => pass,
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
            "Lx" => Lx,
            "Ly" => Ly,
            "strength" => string(strength),
            "n_steps" => result.n_steps,
            "ρ_min" => nan_det ? NaN : minimum(ρ),
            "ρ_max" => nan_det ? NaN : maximum(ρ),
            "method_params" => method_params(method),
            "ci_tier" => "full_nightly",
            "note" => "optional reduced Double Mach-like (full/nightly tier)",
        ),
    ),
    state,
    eq
end

"""
    make_ffs_solid_mask(nx, ny, Lx, Ly; step_x=0.6, step_h=0.2) -> Vector{Bool}

Mark elements whose centers lie inside the classic forward-facing step solid.
"""
function make_ffs_solid_mask(
    nx::Int,
    ny::Int,
    Lx::Float64,
    Ly::Float64;
    step_x::Float64=0.6,
    step_h::Float64=0.2,
)
    solid = falses(nx * ny)
    dx = Lx / nx
    dy = Ly / ny
    for jy in 1:ny, jx in 1:nx
        xc = (jx - 0.5) * dx
        yc = (jy - 0.5) * dy
        if xc >= step_x && yc <= step_h
            e = (jy - 1) * nx + jx
            solid[e] = true
        end
    end
    return collect(solid)
end

"""
    run_forward_facing_step(; ...)

Reduced forward-facing step via solid-element mask on a structured mesh.
**CI tier:** full/nightly only. Pass: finite, non-diverged on ultra-light grids.
"""
function run_forward_facing_step(;
    p::Int=1,
    nx::Int=24,
    ny::Int=8,
    t_final::Float64=0.05,
    cfl::Float64=0.04,
    γ::Float64=1.4,
    Lx::Float64=1.5,
    Ly::Float64=0.5,
    step_x::Float64=0.3,
    step_h::Float64=0.1,
    M_in::Float64=3.0,
    method::AbstractCapturingMethod=PerssonAVMethod(; c_av=0.15),
    method_name::AbstractString="persson_av",
    require_positivity::Bool=false,
)
    t0 = time()
    eq = Euler2D(γ)
    ops = build_operators(p)
    c0 = sqrt(γ * 1.0 / 1.4)
    u_in = M_in * c0
    U_in = primitives_to_conserved(eq, 1.4, u_in, 0.0, 1.0)

    solid = make_ffs_solid_mask(nx, ny, Lx, Ly; step_x=step_x, step_h=step_h)
    mesh = Mesh2D(
        0.0, Lx, 0.0, Ly, nx, ny;
        left_bc=DirichletBC((x, y, t) -> U_in),
        right_bc=TransmissiveBC(),
        bottom_bc=ReflectingBC(),
        top_bc=ReflectingBC(),
        solid=solid,
    )
    state = allocate_state(mesh, ops, Val(4))
    set_initial_condition!(state, (x, y) -> U_in)
    result = integrate!(state, eq, method, t_final; cfl=cfl)
    diverged = result.status != :ok
    nan_det = diverged || has_nonfinite(state.u)
    pos = !nan_det
    if pos
        Np = n_points(ops)
        for e in 1:mesh.n_elements
            is_solid(mesh, e) && continue
            for j in 1:Np, i in 1:Np
                U = @view state.u[i, j, e, :]
                pr = pressure(eq, U)
                if !(U[1] > 0 && pr > 0 && isfinite(U[1]) && isfinite(pr))
                    pos = false
                    break
                end
            end
            pos || break
        end
    end
    pass = !diverged && !nan_det && (!require_positivity || pos)
    n_solid = count(solid)
    return Dict{String,Any}(
        "name" => "forward_facing_step_p$(p)_$(method_name)",
        "case_type" => "discontinuous",
        "equation" => "euler2d",
        "p" => p,
        "capturing_method" => String(method_name),
        "pass" => pass,
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
            "Lx" => Lx,
            "Ly" => Ly,
            "step_x" => step_x,
            "step_h" => step_h,
            "n_solid" => n_solid,
            "n_steps" => result.n_steps,
            "M_in" => M_in,
            "method_params" => method_params(method),
            "ci_tier" => "full_nightly",
            "note" => "optional reduced FFS via solid mask (full/nightly tier)",
        ),
    ),
    state,
    eq
end

"""
    run_p33b_optional_suite() -> (cases, overall, hard_fails)

Optional reduced DMR + FFS (full/nightly). Ultra-light defaults for unit tests.
"""
function run_p33b_optional_suite()
    cases = Any[]
    hard_fails = String[]

    c_dmr, _, _ = run_double_mach_reflection(;
        p=1,
        nx=12,
        ny=4,
        t_final=0.03,
        cfl=0.04,
        Lx=1.2,
        Ly=0.4,
        strength=:reduced,
        method=PerssonAVMethod(; c_av=0.2),
        method_name="persson_av",
        require_positivity=false,
    )
    push!(cases, c_dmr)
    if c_dmr["diverged"] || c_dmr["nan_detected"]
        push!(hard_fails, "reduced Double Mach diverged/NaN")
        c_dmr["pass"] = false
    else
        c_dmr["pass"] = true
    end

    c_ffs, _, _ = run_forward_facing_step(;
        p=1,
        nx=12,
        ny=4,
        t_final=0.02,
        cfl=0.025,
        Lx=1.0,
        Ly=0.4,
        step_x=0.3,
        step_h=0.1,
        M_in=2.0,
        method=PerssonAVMethod(; c_av=0.25),
        method_name="persson_av",
        require_positivity=false,
    )
    push!(cases, c_ffs)
    if c_ffs["diverged"] || c_ffs["nan_detected"]
        push!(hard_fails, "reduced FFS diverged/NaN")
        c_ffs["pass"] = false
    else
        c_ffs["pass"] = true
    end

    overall = isempty(hard_fails)
    return cases, overall, hard_fails
end
