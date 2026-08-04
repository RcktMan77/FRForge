# Optional 2D wall benchmarks (reduced DMR, FFS) + p33b suite.
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
function double_mach_states(eq::Euler2D; strength::Symbol = :reduced)
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
function double_mach_ic(eq::Euler2D, x, y; x_shock0 = 1 / 6, strength::Symbol = :reduced)
    U_post, U_pre = double_mach_states(eq; strength = strength)
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
    p::Int = 1,
    nx::Int = 24,
    ny::Int = 6,
    t_final::Float64 = 0.04,
    cfl::Float64 = 0.04,
    γ::Float64 = 1.4,
    x_shock0::Float64 = 1 / 6,
    Lx::Float64 = 2.0,
    Ly::Float64 = 0.5,
    strength::Symbol = :reduced,
    method::AbstractCapturingMethod = PerssonAVMethod(; c_av = 0.15),
    method_name::AbstractString = "persson_av",
    require_positivity::Bool = false,
    progress_every::Int = 0,
    progress_label::AbstractString = "",
)
    t0 = time()
    eq = Euler2D(γ)
    ops = build_operators(p)
    U_post, U_pre = double_mach_states(eq; strength = strength)
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
        left_bc = left, right_bc = right, bottom_bc = bottom, top_bc = top,
    )
    state = allocate_state(mesh, ops, Val(4))
    set_initial_condition!(
        state, (x, y) -> double_mach_ic(eq, x, y; x_shock0 = x_shock0, strength = strength),
    )
    plabel = isempty(progress_label) ? "double_mach_$(strength)" : progress_label
    result = integrate!(
        state, eq, method, t_final;
        cfl = cfl, progress_every = progress_every, progress_label = plabel,
    )
    diverged = result.status != :ok
    nan_det = diverged || has_nonfinite(state.u)
    pos = !nan_det && positivity_ok(eq, state)
    pass = !diverged && !nan_det && (!require_positivity || pos)
    ρ = @view state.u[:, :, :, 1]
    return case_report_dict(;
    name="double_mach_$(strength)_p$(p)_$(method_name)",
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
    overshoot=0.0,
    metrics=Dict{String,Any}(
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
    step_x::Float64 = 0.6,
    step_h::Float64 = 0.2,
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
    p::Int = 1,
    nx::Int = 24,
    ny::Int = 8,
    t_final::Float64 = 0.05,
    cfl::Float64 = 0.04,
    γ::Float64 = 1.4,
    Lx::Float64 = 1.5,
    Ly::Float64 = 0.5,
    step_x::Float64 = 0.3,
    step_h::Float64 = 0.1,
    M_in::Float64 = 3.0,
    method::AbstractCapturingMethod = PerssonAVMethod(; c_av = 0.15),
    method_name::AbstractString = "persson_av",
    require_positivity::Bool = false,
)
    t0 = time()
    eq = Euler2D(γ)
    ops = build_operators(p)
    c0 = sqrt(γ * 1.0 / 1.4)
    u_in = M_in * c0
    U_in = primitives_to_conserved(eq, 1.4, u_in, 0.0, 1.0)

    solid = make_ffs_solid_mask(nx, ny, Lx, Ly; step_x = step_x, step_h = step_h)
    mesh = Mesh2D(
        0.0, Lx, 0.0, Ly, nx, ny;
        left_bc = DirichletBC((x, y, t) -> U_in),
        right_bc = TransmissiveBC(),
        bottom_bc = ReflectingBC(),
        top_bc = ReflectingBC(),
        solid = solid,
    )
    state = allocate_state(mesh, ops, Val(4))
    set_initial_condition!(state, (x, y) -> U_in)
    result = integrate!(state, eq, method, t_final; cfl = cfl)
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
    return case_report_dict(;
    name="forward_facing_step_p$(p)_$(method_name)",
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
    overshoot=0.0,
    metrics=Dict{String,Any}(
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
        p = 1,
        nx = 12,
        ny = 4,
        t_final = 0.03,
        cfl = 0.04,
        Lx = 1.2,
        Ly = 0.4,
        strength = :reduced,
        method = PerssonAVMethod(; c_av = 0.2),
        method_name = "persson_av",
        require_positivity = false,
    )
    push!(cases, c_dmr)
    if c_dmr["diverged"] || c_dmr["nan_detected"]
        push!(hard_fails, "reduced Double Mach diverged/NaN")
        c_dmr["pass"] = false
    else
        c_dmr["pass"] = true
    end

    c_ffs, _, _ = run_forward_facing_step(;
        p = 1,
        nx = 12,
        ny = 4,
        t_final = 0.02,
        cfl = 0.025,
        Lx = 1.0,
        Ly = 0.4,
        step_x = 0.3,
        step_h = 0.1,
        M_in = 2.0,
        method = PerssonAVMethod(; c_av = 0.25),
        method_name = "persson_av",
        require_positivity = false,
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
