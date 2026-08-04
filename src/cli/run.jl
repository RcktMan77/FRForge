# `frforge run` — single-case driver.
#
# Each CLI case is a small helper (`_run_case_*`) so the dispatcher stays thin.
# Numerics match the previous monolithic `cli_run` (same meshes, CFLs, defaults).

# ---------------------------------------------------------------------------
# Mesh / state setup (CLI only; behavior-preserving patterns)
# ---------------------------------------------------------------------------

function _cli_mesh1d_periodic(xL::Real, xR::Real, Ne::Int)
    return Mesh1D(Float64(xL), Float64(xR), Ne; left_bc=PeriodicBC(), right_bc=PeriodicBC())
end

function _cli_mesh1d_transmissive(xL::Real, xR::Real, Ne::Int)
    return Mesh1D(
        Float64(xL),
        Float64(xR),
        Ne;
        left_bc=TransmissiveBC(),
        right_bc=TransmissiveBC(),
    )
end

function _cli_mesh2d_unit(Ne::Int; transmissive::Bool=false)
    if transmissive
        return Mesh2D(
            0.0,
            1.0,
            0.0,
            1.0,
            Ne,
            Ne;
            left_bc=TransmissiveBC(),
            right_bc=TransmissiveBC(),
            bottom_bc=TransmissiveBC(),
            top_bc=TransmissiveBC(),
        )
    end
    return Mesh2D(0.0, 1.0, 0.0, 1.0, Ne, Ne)
end

"""When CLI default t_final is 1.0, use case-specific default; else honor the user value."""
function _cli_t_final(opts::AbstractDict, case_default::Float64)
    return opts["t_final"] == 1.0 ? case_default : opts["t_final"]
end

function _cli_state1d(mesh, ops, scheme, ::Val{Neq}, u0) where {Neq}
    state = allocate_state(mesh, ops, Val(Neq); scheme=scheme)
    set_initial_condition!(state, u0)
    return state
end

function _cli_state2d(mesh, ops, scheme, ::Val{Neq}, u0) where {Neq}
    state = allocate_state(mesh, ops, Val(Neq); scheme=scheme)
    set_initial_condition!(state, u0)
    return state
end

"""Optionally write high-order VTU if --output is set."""
function _maybe_write_vtu(opts, state, eq, method=nothing)
    out = opts["output"]
    if !isempty(out)
        if get(opts, "vtk_diagnostics", false) &&
           state isa SolutionState2D &&
           method !== nothing
            write_vtu_high_order_with_capturing(out, state, eq, method)
            println(
                "VTU written (with sensor/av diagnostics): $(abspath(out))  (ParaView ≥ 5.5)",
            )
        else
            write_vtu_high_order(out, state, eq)
            println("VTU written: $(abspath(out))  (open in ParaView ≥ 5.5)")
        end
        return true
    end
    return false
end

# ---------------------------------------------------------------------------
# Per-case helpers (return process exit code 0/1)
# ---------------------------------------------------------------------------

function _run_case_advection_sine(opts, ops, scheme, method)
    p, Ne, t_final, cfl = opts["p"], opts["ne"], opts["t_final"], opts["cfl"]
    a = opts["a"]
    eq = LinearAdvection1D(a)
    mesh = _cli_mesh1d_periodic(0.0, 1.0, Ne)
    state = _cli_state1d(mesh, ops, scheme, Val(1), x -> sin(2π * x))
    M0 = discrete_mass(state, 1)
    result = integrate!(state, eq, method, t_final; cfl=cfl)
    MT = discrete_mass(state, 1)
    err = l2_error(state, x -> sin(2π * (x - a * t_final)), 1)
    println(
        "case=advection_sine p=$p ne=$Ne a=$a t_final=$t_final cfl=$cfl method=$(opts["method"]) scheme=$(scheme_dict(scheme))",
    )
    println("status=$(result.status) n_steps=$(result.n_steps) t=$(result.t)")
    println("L2_error=$err mass_change=$(MT - M0)")
    _maybe_write_vtu(opts, state, eq)
    return result.status == :ok ? 0 : 1
end

function _run_case_burgers_square(opts, ops, scheme, method)
    p, Ne, t_final, cfl = opts["p"], opts["ne"], opts["t_final"], opts["cfl"]
    eq = Burgers1D()
    mesh = _cli_mesh1d_periodic(0.0, 1.0, Ne)
    state = _cli_state1d(mesh, ops, scheme, Val(1), x -> burgers_square_ic(x))
    u0_min, u0_max = solution_extrema(state, 1)
    M0 = discrete_mass(state, 1)
    result = integrate!(state, eq, method, t_final; cfl=cfl)
    MT = discrete_mass(state, 1)
    u_min, u_max = solution_extrema(state, 1)
    _, _, η = overshoot_metric(u_min, u_max, u0_min, u0_max)
    println(
        "case=burgers_square p=$p ne=$Ne t_final=$t_final cfl=$cfl method=$(opts["method"]) scheme=$(scheme_dict(scheme))",
    )
    println("status=$(result.status) n_steps=$(result.n_steps) t=$(result.t)")
    println("mass_change=$(MT - M0) u_range=[$u_min, $u_max] overshoot=$η")
    if opts["method"] == "null"
        println("note: pure high-order FR is expected to oscillate near the discontinuity")
    end
    _maybe_write_vtu(opts, state, eq)
    return result.status == :ok ? 0 : 1
end

function _run_case_euler_density_wave(opts, ops, scheme, method)
    p, Ne, t_final, cfl = opts["p"], opts["ne"], opts["t_final"], opts["cfl"]
    eq = Euler1D(1.4)
    mesh = _cli_mesh1d_periodic(0.0, 1.0, Ne)
    state = _cli_state1d(mesh, ops, scheme, Val(3), x -> euler_density_wave_conserved(eq, x, 0.0))
    M0 = discrete_mass(state, 1)
    result = integrate!(state, eq, method, t_final; cfl=cfl)
    MT = discrete_mass(state, 1)
    err = l2_error(state, x -> euler_density_wave_conserved(eq, x, t_final), 1)
    println(
        "case=euler_density_wave p=$p ne=$Ne t_final=$t_final cfl=$cfl method=$(opts["method"]) scheme=$(scheme_dict(scheme))",
    )
    println("status=$(result.status) n_steps=$(result.n_steps) t=$(result.t)")
    println("L2_density=$err mass_change=$(MT - M0) positivity=$(positivity_ok(eq, state))")
    _maybe_write_vtu(opts, state, eq)
    return result.status == :ok ? 0 : 1
end

function _run_case_sod(opts, ops, scheme, method)
    p, Ne, cfl = opts["p"], opts["ne"], opts["cfl"]
    tf = _cli_t_final(opts, 0.2)
    # Run solver again for VTU (run_sod returns metrics only)
    eq = Euler1D(1.4)
    mesh = _cli_mesh1d_transmissive(0.0, 1.0, Ne)
    state = _cli_state1d(mesh, ops, scheme, Val(3), x -> sod_ic(eq, x))
    result = integrate!(state, eq, method, tf; cfl=min(cfl, 0.15))
    c = run_sod(;
        p=p,
        n_elements=Ne,
        t_final=tf,
        cfl=cfl,
        method=method,
        method_name=opts["method"],
    )
    println(
        "case=sod p=$p ne=$Ne t_final=$tf method=$(opts["method"]) scheme=$(scheme_dict(scheme))",
    )
    println(
        "status=$(c["pass"] ? "ok" : "fail") η=$(c["overshoot"]) δ=$(c["shock_thickness"]) Dex=$(c["excess_dissipation"])",
    )
    println("L1_vs_exact=$(c["l1_error_vs_reference"]) positivity=$(c["positivity_ok"])")
    if result.status == :ok
        _maybe_write_vtu(opts, state, eq)
    end
    return c["pass"] ? 0 : 1
end

function _run_case_shu_osher(opts, ops, scheme, method)
    p, Ne, cfl = opts["p"], opts["ne"], opts["cfl"]
    tf = _cli_t_final(opts, 1.8)
    eq = Euler1D(1.4)
    mesh = _cli_mesh1d_transmissive(0.0, 10.0, Ne)
    state = _cli_state1d(mesh, ops, scheme, Val(3), x -> shu_osher_ic(eq, x))
    result = integrate!(state, eq, method, tf; cfl=cfl)
    c = run_shu_osher(;
        p=p,
        n_elements=Ne,
        t_final=tf,
        cfl=cfl,
        method=method,
        method_name=opts["method"],
    )
    println(
        "case=shu_osher p=$p ne=$Ne t_final=$tf method=$(opts["method"]) scheme=$(scheme_dict(scheme))",
    )
    println(
        "status=$(c["pass"] ? "ok" : "fail") η=$(c["overshoot"]) δ=$(c["shock_thickness"]) Dex=$(c["excess_dissipation"])",
    )
    println("positivity=$(c["positivity_ok"])")
    if result.status == :ok
        _maybe_write_vtu(opts, state, eq)
    end
    return c["pass"] ? 0 : 1
end

function _run_case_advection2d(opts, ops, scheme, method)
    p, Ne, t_final, cfl = opts["p"], opts["ne"], opts["t_final"], opts["cfl"]
    ax, ay = 1.0, 0.5
    eq = LinearAdvection2D(ax, ay)
    mesh = _cli_mesh2d_unit(Ne)
    state = _cli_state2d(mesh, ops, scheme, Val(1), (x, y) -> sin(2π * x) * sin(2π * y))
    result = integrate!(state, eq, method, t_final; cfl=cfl)
    err = l2_error(
        state,
        (x, y) -> sin(2π * (x - ax * t_final)) * sin(2π * (y - ay * t_final)),
        1,
    )
    println(
        "case=advection2d p=$p ne=$(Ne)x$(Ne) t_final=$t_final L2=$err status=$(result.status) scheme=$(scheme_dict(scheme))",
    )
    _maybe_write_vtu(opts, state, eq)
    return result.status == :ok ? 0 : 1
end

function _run_case_euler2d_wave(opts, ops, scheme, method)
    p, Ne, t_final, cfl = opts["p"], opts["ne"], opts["t_final"], opts["cfl"]
    eq = Euler2D(1.4)
    mesh = _cli_mesh2d_unit(Ne)
    state = _cli_state2d(
        mesh,
        ops,
        scheme,
        Val(4),
        (x, y) ->
            primitives_to_conserved(
                eq,
                1.0 + 0.2 * sin(2π * x) * sin(2π * y),
                1.0,
                1.0,
                1.0,
            ),
    )
    result = integrate!(state, eq, method, t_final; cfl=min(cfl, 0.15))
    println(
        "case=euler2d_wave p=$p ne=$(Ne)x$(Ne) t_final=$t_final status=$(result.status) pos=$(positivity_ok(eq, state)) scheme=$(scheme_dict(scheme))",
    )
    _maybe_write_vtu(opts, state, eq, method)
    return result.status == :ok ? 0 : 1
end

function _run_case_euler2d_jump(opts, ops, scheme, method)
    p, Ne, cfl = opts["p"], opts["ne"], opts["cfl"]
    tf = _cli_t_final(opts, 0.05)
    eq = Euler2D(1.4)
    mesh = _cli_mesh2d_unit(Ne; transmissive=true)
    state = _cli_state2d(
        mesh,
        ops,
        scheme,
        Val(4),
        (x, y) ->
            x < 0.5 ? primitives_to_conserved(eq, 1.0, 0.0, 0.0, 1.0) :
            primitives_to_conserved(eq, 0.125, 0.0, 0.0, 0.1),
    )
    result = integrate!(state, eq, method, tf; cfl=min(cfl, 0.1))
    println(
        "case=euler2d_jump p=$p ne=$(Ne)x$(Ne) t_final=$tf status=$(result.status) pos=$(positivity_ok(eq, state)) scheme=$(scheme_dict(scheme))",
    )
    _maybe_write_vtu(opts, state, eq, method)
    return result.status == :ok ? 0 : 1
end

function _run_case_riemann2d(opts, ops, scheme, method)
    # ops/scheme unused: run_euler2d_riemann builds its own operators (same as before)
    p, Ne, cfl = opts["p"], opts["ne"], opts["cfl"]
    tf = _cli_t_final(opts, 0.12)
    ny = opts["ny"] > 0 ? opts["ny"] : Ne
    c, state, eq = run_euler2d_riemann(;
        p=p,
        nx=Ne,
        ny=ny,
        t_final=tf,
        cfl=min(cfl, 0.08),
        config=:cfg6,
        method=method,
        method_name=opts["method"],
    )
    println(
        "case=riemann2d cfg6 p=$p ne=$(Ne)x$(ny) t_final=$tf method=$(opts["method"]) pass=$(c["pass"]) pos=$(c["positivity_ok"])",
    )
    if !c["diverged"]
        _maybe_write_vtu(opts, state, eq, method)
    end
    return c["pass"] ? 0 : 1
end

function _run_case_double_mach(opts, ops, scheme, method)
    p, Ne, cfl = opts["p"], opts["ne"], opts["cfl"]
    tf = _cli_t_final(opts, 0.08)
    ny = opts["ny"] > 0 ? opts["ny"] : max(div(Ne, 3), 4)
    c, state, eq = run_double_mach_reflection(;
        p=p,
        nx=Ne,
        ny=ny,
        t_final=tf,
        cfl=min(cfl, 0.05),
        Lx=1.5,
        Ly=0.5,
        strength=:reduced,
        method=method,
        method_name=opts["method"],
        require_positivity=false,
    )
    println(
        "case=double_mach reduced p=$p ne=$(Ne)x$(ny) t_final=$tf method=$(opts["method"]) pass=$(c["pass"]) pos=$(c["positivity_ok"])",
    )
    if !c["diverged"]
        _maybe_write_vtu(opts, state, eq, method)
    end
    return c["pass"] ? 0 : 1
end

# ---------------------------------------------------------------------------
# Arg parse + dispatcher
# ---------------------------------------------------------------------------

function _parse_run_args(args)
    s = ArgParseSettings(
        description="Run a single FRForge case.",
        prog="frforge run",
    )
    @add_arg_table! s begin
        "--case"
        help = CLI_RUN_CASE_HELP
        default = "advection_sine"
        "--p"
        help = "Polynomial degree"
        arg_type = Int
        default = 3
        "--ne"
        help = "Number of elements"
        arg_type = Int
        default = 16
        "--t-final"
        help = "Final time"
        arg_type = Float64
        default = 1.0
        dest_name = "t_final"
        "--cfl"
        help = "CFL number"
        arg_type = Float64
        default = 0.2
        "--a"
        help = "Advection speed"
        arg_type = Float64
        default = 1.0
        "--method"
        help = "Capturing method name (see frforge --help for registry)"
        default = "null"
        "--points"
        help = "Solution points: gl (default) | gll"
        default = "gl"
        "--flux"
        help = "Numerical flux: rusanov (default) | hllc"
        default = "rusanov"
        "--time"
        help = "Time integrator: ssp_rk3 (default) | ssp_rk2"
        default = "ssp_rk3"
        "--output"
        help = "Optional high-order VTU output path (ParaView)"
        default = ""
        dest_name = "output"
        "--vtk-diagnostics"
        help = "With --output on 2D Euler: write sensor + av cell/point fields (docs only)"
        action = :store_true
        dest_name = "vtk_diagnostics"
        "--ny"
        help = "Y elements for 2D non-square cases (default: --ne)"
        arg_type = Int
        default = 0
    end
    return parse_args(args, s)
end

const _CLI_RUN_DISPATCH = Dict{String,Function}(
    "advection_sine" => _run_case_advection_sine,
    "burgers_square" => _run_case_burgers_square,
    "euler_density_wave" => _run_case_euler_density_wave,
    "sod" => _run_case_sod,
    "shu_osher" => _run_case_shu_osher,
    "advection2d" => _run_case_advection2d,
    "euler2d_wave" => _run_case_euler2d_wave,
    "euler2d_jump" => _run_case_euler2d_jump,
    "riemann2d" => _run_case_riemann2d,
    "double_mach" => _run_case_double_mach,
)

function cli_run(args)
    opts = _parse_run_args(args)
    case = opts["case"]
    runner = get(_CLI_RUN_DISPATCH, case, nothing)
    if runner === nothing
        println(stderr, "Unknown case: $case")
        println(stderr, "  ", CLI_RUN_CASE_HELP)
        return 2
    end
    method = get_capturing_method(opts["method"])
    scheme = _scheme_from_opts(opts)
    ops = build_operators(opts["p"]; points=scheme.points)
    return runner(opts, ops, scheme, method)
end
