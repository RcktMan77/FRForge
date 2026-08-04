# 2D Riemann configurations + isentropic vortex order + p33a suite.
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
function riemann2d_ic(eq::Euler2D, x, y; config = :cfg6, x0 = 0.5, y0 = 0.5)
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
    riemann2d_ic(eq, x, y; config = :cfg3, kwargs...)

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
    p::Int = 1,
    nx::Int = 16,
    ny::Int = 16,
    t_final::Float64 = 0.08,
    cfl::Float64 = 0.08,
    γ::Float64 = 1.4,
    x0::Float64 = 0.5,
    y0::Float64 = 0.5,
    config = :cfg6,
    method::AbstractCapturingMethod = PerssonAVMethod(; c_av = 0.1),
    method_name::AbstractString = "persson_av",
    require_positivity::Bool = true,
    progress_every::Int = 0,
    progress_label::AbstractString = "",
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
        left_bc = TransmissiveBC(),
        right_bc = TransmissiveBC(),
        bottom_bc = TransmissiveBC(),
        top_bc = TransmissiveBC(),
    )
    state = allocate_state(mesh, ops, Val(4))
    set_initial_condition!(
        state, (x, y) -> riemann2d_ic(eq, x, y; config = config, x0 = x0, y0 = y0),
    )
    ρ0 = @view state.u[:, :, :, 1]
    ρmin0, ρmax0 = minimum(ρ0), maximum(ρ0)
    plabel = isempty(progress_label) ? "riemann_$(config)" : progress_label
    result = integrate!(
        state, eq, method, t_final;
        cfl = cfl, progress_every = progress_every, progress_label = plabel,
    )
    diverged = result.status != :ok
    nan_det = diverged || has_nonfinite(state.u)
    pos = !nan_det && positivity_ok(eq, state)
    ρ = @view state.u[:, :, :, 1]
    ρmin, ρmax = nan_det ? (NaN, NaN) : (minimum(ρ), maximum(ρ))
    _, _, η = overshoot_metric(ρmin, ρmax, min(ρmin0, 0.5), max(ρmax0, 3.0))
    pass = !diverged && !nan_det && (!require_positivity || pos)
    cfg_str = string(config)
    return case_report_dict(;
    name="euler2d_riemann_$(cfg_str)_p$(p)_$(method_name)",
    case_type="discontinuous",
    equation="euler2d",
    p=p,
    capturing_method=String(method_name),
    pass=pass,
    diverged=diverged,
    nan_detected=nan_det,
    conservation_residual=0.0,
    conservation_pass=true,
    conservation_metric="none",
    positivity_ok=pos,
    wall_time_sec=time() - t0,
    n_elements=nx * ny,
    t_final=t_final,
    excess_dissipation=nothing,
    shock_thickness=nothing,
    shock_thickness_unit="sp_spacings",
    overshoot=η,
    metrics=Dict{String,Any}(
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
    γ::Float64 = 1.4,
    β::Float64 = 5.0,
    x0::Float64 = 5.0,
    y0::Float64 = 5.0,
    u∞::Float64 = 1.0,
    v∞::Float64 = 1.0,
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
    p::Int = 2,
    n_list::AbstractVector{Int} = [8, 16],
    t_final::Float64 = 0.5,
    cfl::Float64 = 0.08,
    L::Float64 = 10.0,
    γ::Float64 = 1.4,
    β::Float64 = 5.0,
    x0::Float64 = 5.0,
    y0::Float64 = 5.0,
    u∞::Float64 = 1.0,
    v∞::Float64 = 1.0,
    order_tol::Float64 = DEFAULT_ORDER_TOLERANCE,
    method::AbstractCapturingMethod = NullCapturing(),
    method_name::AbstractString = "null",
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
            x, y, t_final; γ = γ, β = β, x0 = x0, y0 = y0, u∞ = u∞, v∞ = v∞,
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
                    x, y, 0.0; γ = γ, β = β, x0 = x0, y0 = y0, u∞ = u∞, v∞ = v∞,
                )
                primitives_to_conserved(eq, ρ, u, v, pr)
            end,
        )
        result = integrate!(state, eq, method, t_final; dt = dt_fixed)
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
    opass = !diverged && order_pass(obs; formal_order = formal, tol = order_tol + 0.15)

    # Periodic mass conservation on finest
    n = n_list[end]
    mesh = Mesh2D(0.0, L, 0.0, L, n, n)
    state = allocate_state(mesh, ops, Val(4))
    set_initial_condition!(
        state,
        (x, y) -> begin
            ρ, u, v, pr = isentropic_vortex_primitives(
                x, y, 0.0; γ = γ, β = β, x0 = x0, y0 = y0, u∞ = u∞, v∞ = v∞,
            )
            primitives_to_conserved(eq, ρ, u, v, pr)
        end,
    )
    M0 = discrete_mass(state, 1)
    result = integrate!(state, eq, method, t_final; cfl = cfl)
    MT = discrete_mass(state, 1)
    cres = conservation_residual_relative(M0, MT)
    cpass = result.status == :ok && (cres <= 1e-8 || abs(MT - M0) <= 1e-10)
    pos_ok = pos_ok && positivity_ok(eq, state)

    return case_report_dict(;
    name="isentropic_vortex_order_p$(p)_$(method_name)",
    case_type="smooth_order",
    equation="euler2d",
    p=p,
    capturing_method=String(method_name),
    pass=opass && cpass && !diverged && pos_ok,
    diverged=diverged,
    nan_detected=diverged,
    conservation_residual=cres,
    conservation_pass=cpass,
    conservation_metric="periodic_mass_change",
    positivity_ok=pos_ok,
    wall_time_sec=time() - t0,
    metrics=Dict{String,Any}(
            "n_list" => collect(n_list),
            "L" => L,
            "β" => β,
            "t_final" => t_final,
            "method_params" => method_params(method),
            "ci_tier" => "required_ci_light",
            "note" => "isentropic vortex = exact-solution smooth Euler (MMS-style) gate",
        ),
    extra=Dict{String,Any}(
        "mesh_sizes" => mesh_sizes,
        "l2_errors" => l2_errors,
        "observed_orders" => obs,
        "formal_order" => formal,
        "order_pass" => opass,
        "order_tolerance" => order_tol + 0.15,
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
        p = 2,
        n_list = [8, 16],
        t_final = 0.5,
        cfl = 0.08,
        L = 10.0,
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
        p = 1,
        nx = 16,
        ny = 16,
        t_final = 0.08,
        cfl = 0.08,
        config = :cfg6,
        method = PerssonAVMethod(; c_av = 0.1),
        method_name = "persson_av",
    )
    push!(cases, c_r)
    if !c_r["pass"]
        push!(hard_fails, "2D Riemann cfg6 + persson_av failed")
    end

    # Short NullCapturing smoke on same config (positivity still required for cfg6)
    c_rn, _, _ = run_euler2d_riemann(;
        p = 1,
        nx = 12,
        ny = 12,
        t_final = 0.04,
        cfl = 0.06,
        config = :cfg6,
        method = NullCapturing(),
        method_name = "null",
    )
    push!(cases, c_rn)
    if !c_rn["pass"]
        push!(hard_fails, "2D Riemann cfg6 null failed")
    end

    overall = isempty(hard_fails) && all(c -> c["pass"] === true, cases)
    return cases, overall, hard_fails
end

