# CLI entry points for the `frforge` binary.

function _print_usage(io=stderr)
    println(io, "FRForge — high-order Flux Reconstruction laboratory")
    println(io, "Usage: frforge {test|run|invent|score} [options]")
    println(io)
    println(io, "Commands:")
    println(io, "  test    Run verification suite and emit JSON report")
    println(io, "  run     Run a single case")
    println(io, "  invent  Run quant suite for a method vs baseline and classify candidate")
    println(io, "  score   Classify two existing JSON reports (method vs baseline)")
    println(io)
    println(io, "Examples:")
    println(io, "  frforge test --suite quant --method persson_av --report results/m5/report.json")
    println(io, "  frforge invent --method scaled_persson --baseline persson_av")
    println(io, "  frforge score --method-report a.json --baseline-report b.json")
    println(io, "  frforge run --case sod --p 2 --ne 64 --method persson_av")
    println(io, describe_methods())
end

function _parse_test_args(args)
    s = ArgParseSettings(
        description="Run FRForge verification and write a JSON report.",
        prog="frforge test",
    )
    @add_arg_table! s begin
        "--report", "-r"
        help = "Path for the JSON report"
        default = "results/report.json"
        dest_name = "report"
        "--suite"
        help = "Suite: smoke|advection|burgers|euler|capturing|quant|m5|full"
        default = "advection"
        "--method"
        help = "Capturing method name (see frforge --help for registry)"
        default = "null"
    end
    return parse_args(args, s)
end

function _parse_invent_args(args)
    s = ArgParseSettings(
        description="Invent: run quant suite for method vs baseline and classify.",
        prog="frforge invent",
    )
    @add_arg_table! s begin
        "--method", "-m"
        help = "Candidate method name (registered)"
        required = true
        "--baseline", "-b"
        help = "Baseline method name"
        default = "persson_av"
        "--report-dir"
        help = "Directory for invent JSON artifacts"
        default = "results/invent"
        dest_name = "report_dir"
        "--delta"
        help = "Composite margin δ_score for promising status"
        arg_type = Float64
        default = DEFAULT_SCORE_MARGIN
        "--vtk-produced"
        help = "Set true if HO VTK was produced (accepted_candidate)"
        action = :store_true
        dest_name = "vtk_produced"
    end
    return parse_args(args, s)
end

function _parse_score_args(args)
    s = ArgParseSettings(
        description="Score two existing reports (method vs baseline).",
        prog="frforge score",
    )
    @add_arg_table! s begin
        "--method-report"
        help = "Path to method report JSON"
        required = true
        dest_name = "method_report"
        "--baseline-report"
        help = "Path to baseline report JSON"
        required = true
        dest_name = "baseline_report"
        "--delta"
        help = "Composite margin δ_score"
        arg_type = Float64
        default = DEFAULT_SCORE_MARGIN
        "--vtk-produced"
        action = :store_true
        dest_name = "vtk_produced"
        "--output", "-o"
        help = "Optional path for comparison JSON"
        default = ""
    end
    return parse_args(args, s)
end

function _parse_run_args(args)
    s = ArgParseSettings(
        description="Run a single FRForge case.",
        prog="frforge run",
    )
    @add_arg_table! s begin
        "--case"
        help = "Case: advection_sine|…|sod|shu_osher|advection2d|euler2d_wave|euler2d_jump"
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
        help = "Capturing method: null | persson_av | scaled_persson"
        default = "null"
        "--output"
        help = "Optional high-order VTU output path (ParaView)"
        default = ""
        dest_name = "output"
    end
    return parse_args(args, s)
end

"""
    cli_test(opts) -> Int

Run verification suite and write schema v1 JSON report.
Suites: smoke | advection/m1 | burgers/m2 | euler/m3 | full.
"""
function cli_test(opts::AbstractDict)
    t0 = time()
    report_path = opts["report"]
    suite = opts["suite"]
    method = opts["method"]

    cases = Any[]
    overall_pass = true
    diverged = false
    nan_detected = false
    hard_fails = String[]

    if suite == "smoke"
        # M0 skeleton
        cases = Any[]
        overall_pass = true
    elseif suite in ("advection", "m1")
        cases, overall_pass, hard_fails = run_m1_advection_suite()
        diverged = any(c -> get(c, "diverged", false) === true, cases)
        nan_detected = any(c -> get(c, "nan_detected", false) === true, cases)
    elseif suite in ("burgers", "m2")
        cases, overall_pass, hard_fails = run_m2_burgers_suite()
        diverged = any(c -> get(c, "diverged", false) === true, cases)
        nan_detected = any(c -> get(c, "nan_detected", false) === true, cases)
    elseif suite in ("euler", "m3")
        cases, overall_pass, hard_fails = run_m3_euler_suite()
        diverged = any(c -> get(c, "diverged", false) === true, cases)
        nan_detected = any(c -> get(c, "nan_detected", false) === true, cases)
    elseif suite in ("capturing", "m4")
        cases, overall_pass, hard_fails = run_m4_capturing_suite()
        diverged = any(c -> get(c, "diverged", false) === true, cases)
        nan_detected = any(c -> get(c, "nan_detected", false) === true, cases)
        method = "persson_av"
    elseif suite in ("quant", "m5")
        cases, overall_pass, hard_fails = run_m5_quant_suite(; method_name=method)
        diverged = any(c -> get(c, "diverged", false) === true, cases)
        nan_detected = any(c -> get(c, "nan_detected", false) === true, cases)
    elseif suite in ("2d", "m8")
        cases, overall_pass, hard_fails, _, _ = run_m8_2d_suite()
        diverged = any(c -> get(c, "diverged", false) === true, cases)
        nan_detected = any(c -> get(c, "nan_detected", false) === true, cases)
    elseif suite == "full"
        c1, p1, f1 = run_m1_advection_suite()
        c2, p2, f2 = run_m2_burgers_suite()
        c3, p3, f3 = run_m3_euler_suite()
        c4, p4, f4 = run_m4_capturing_suite()
        c5, p5, f5 = run_m5_quant_suite(; method_name=method)
        c8, p8, f8, _, _ = run_m8_2d_suite()
        cases = vcat(c1, c2, c3, c4, c5, c8)
        hard_fails = vcat(f1, f2, f3, f4, f5, f8)
        overall_pass = p1 && p2 && p3 && p4 && p5 && p8
        diverged = any(c -> get(c, "diverged", false) === true, cases)
        nan_detected = any(c -> get(c, "nan_detected", false) === true, cases)
    else
        println(stderr, "Unknown suite: $suite (use smoke|…|quant|2d|full)")
        return 2
    end

    mparams = Dict{String,Any}()
    try
        mparams = method_params(get_capturing_method(method))
    catch
        mparams = Dict{String,Any}("method" => method)
    end

    report = write_report_skeleton(
        report_path;
        command="test",
        suite=suite,
        method_name=method,
        method_params=mparams,
        baseline_name=nothing,
        overall_pass=overall_pass,
        diverged=diverged,
        nan_detected=nan_detected,
        wall_time_sec=time() - t0,
        hard_gate_failures=hard_fails,
        cases=cases,
    )

    errs = validate_report_keys(report)
    if !isempty(errs)
        println(stderr, "Report validation failed:")
        for e in errs
            println(stderr, "  - ", e)
        end
        return 1
    end

    println("FRForge test suite=$(suite) method=$(method)")
    println("Report written: $(abspath(report_path))")
    sc = report["summary"]["scores"]
    println(
        "schema_version=$(report["schema_version"]) overall_pass=$(report["overall_pass"]) n_cases=$(report["summary"]["n_cases"])",
    )
    if sc["composite"] !== nothing
        println(
            "scores: order=$(sc["order_preservation"]) dissip=$(sc["dissipation"]) shock=$(sc["shock_quality"]) robust=$(sc["robustness"]) composite=$(sc["composite"])",
        )
    end
    if !isempty(hard_fails)
        println(stderr, "Hard gate failures:")
        for f in hard_fails
            println(stderr, "  - ", f)
        end
    end
    for c in cases
        status = c["pass"] ? "PASS" : "FAIL"
        extra = ""
        if c["case_type"] == "smooth_order"
            extra = " orders=$(c["observed_orders"]) formal=$(c["formal_order"])"
        elseif c["case_type"] == "discontinuous" && haskey(c, "overshoot")
            dex = get(c, "excess_dissipation", nothing)
            δ = get(c, "shock_thickness", nothing)
            extra = " η=$(c["overshoot"]) δ=$(δ) Dex=$(dex)"
        elseif haskey(c, "conservation_residual")
            extra = " cons_res=$(c["conservation_residual"])"
        end
        println("  [$status] $(c["name"])$extra")
    end

    return overall_pass ? 0 : 1
end

"""Optionally write high-order VTU if --output is set."""
function _maybe_write_vtu(opts, state, eq)
    out = opts["output"]
    if !isempty(out)
        write_vtu_high_order(out, state, eq)
        println("VTU written: $(abspath(out))  (open in ParaView ≥ 5.5)")
        return true
    end
    return false
end

function cli_run(args)
    opts = _parse_run_args(args)
    case = opts["case"]
    p = opts["p"]
    Ne = opts["ne"]
    t_final = opts["t_final"]
    cfl = opts["cfl"]
    method = get_capturing_method(opts["method"])

    if case == "advection_sine"
        a = opts["a"]
        eq = LinearAdvection1D(a)
        ops = build_operators(p)
        mesh = Mesh1D(0.0, 1.0, Ne; left_bc=PeriodicBC(), right_bc=PeriodicBC())
        state = allocate_state(mesh, ops, Val(1))
        set_initial_condition!(state, x -> sin(2π * x))
        M0 = discrete_mass(state, 1)
        result = ssp_rk3!(state, eq, method, t_final; cfl=cfl)
        MT = discrete_mass(state, 1)
        err = l2_error(state, x -> sin(2π * (x - a * t_final)), 1)
        println("case=advection_sine p=$p ne=$Ne a=$a t_final=$t_final cfl=$cfl method=$(opts["method"])")
        println("status=$(result.status) n_steps=$(result.n_steps) t=$(result.t)")
        println("L2_error=$err mass_change=$(MT - M0)")
        _maybe_write_vtu(opts, state, eq)
        return result.status == :ok ? 0 : 1
    elseif case == "burgers_square"
        eq = Burgers1D()
        ops = build_operators(p)
        mesh = Mesh1D(0.0, 1.0, Ne; left_bc=PeriodicBC(), right_bc=PeriodicBC())
        state = allocate_state(mesh, ops, Val(1))
        set_initial_condition!(state, x -> burgers_square_ic(x))
        u0_min, u0_max = solution_extrema(state, 1)
        M0 = discrete_mass(state, 1)
        result = ssp_rk3!(state, eq, method, t_final; cfl=cfl)
        MT = discrete_mass(state, 1)
        u_min, u_max = solution_extrema(state, 1)
        _, _, η = overshoot_metric(u_min, u_max, u0_min, u0_max)
        println("case=burgers_square p=$p ne=$Ne t_final=$t_final cfl=$cfl method=$(opts["method"])")
        println("status=$(result.status) n_steps=$(result.n_steps) t=$(result.t)")
        println("mass_change=$(MT - M0) u_range=[$u_min, $u_max] overshoot=$η")
        if opts["method"] == "null"
            println("note: pure high-order FR is expected to oscillate near the discontinuity")
        end
        _maybe_write_vtu(opts, state, eq)
        return result.status == :ok ? 0 : 1
    elseif case == "euler_density_wave"
        eq = Euler1D(1.4)
        ops = build_operators(p)
        mesh = Mesh1D(0.0, 1.0, Ne; left_bc=PeriodicBC(), right_bc=PeriodicBC())
        state = allocate_state(mesh, ops, Val(3))
        set_initial_condition!(state, x -> euler_density_wave_conserved(eq, x, 0.0))
        M0 = discrete_mass(state, 1)
        result = ssp_rk3!(state, eq, method, t_final; cfl=cfl)
        MT = discrete_mass(state, 1)
        err = l2_error(state, x -> euler_density_wave_conserved(eq, x, t_final), 1)
        println("case=euler_density_wave p=$p ne=$Ne t_final=$t_final cfl=$cfl method=$(opts["method"])")
        println("status=$(result.status) n_steps=$(result.n_steps) t=$(result.t)")
        println("L2_density=$err mass_change=$(MT - M0) positivity=$(positivity_ok(eq, state))")
        _maybe_write_vtu(opts, state, eq)
        return result.status == :ok ? 0 : 1
    elseif case == "sod"
        tf = opts["t_final"] == 1.0 ? 0.2 : opts["t_final"]
        # Run solver again for VTU (run_sod returns metrics only)
        eq = Euler1D(1.4)
        ops = build_operators(p)
        mesh = Mesh1D(0.0, 1.0, Ne; left_bc=TransmissiveBC(), right_bc=TransmissiveBC())
        state = allocate_state(mesh, ops, Val(3))
        set_initial_condition!(state, x -> sod_ic(eq, x))
        result = ssp_rk3!(state, eq, method, tf; cfl=min(cfl, 0.15))
        c = run_sod(;
            p=p,
            n_elements=Ne,
            t_final=tf,
            cfl=cfl,
            method=method,
            method_name=opts["method"],
        )
        println("case=sod p=$p ne=$Ne t_final=$tf method=$(opts["method"])")
        println("status=$(c["pass"] ? "ok" : "fail") η=$(c["overshoot"]) δ=$(c["shock_thickness"]) Dex=$(c["excess_dissipation"])")
        println("L1_vs_exact=$(c["l1_error_vs_reference"]) positivity=$(c["positivity_ok"])")
        if result.status == :ok
            _maybe_write_vtu(opts, state, eq)
        end
        return c["pass"] ? 0 : 1
    elseif case == "shu_osher"
        tf = opts["t_final"] == 1.0 ? 1.8 : opts["t_final"]
        eq = Euler1D(1.4)
        ops = build_operators(p)
        mesh = Mesh1D(0.0, 10.0, Ne; left_bc=TransmissiveBC(), right_bc=TransmissiveBC())
        state = allocate_state(mesh, ops, Val(3))
        set_initial_condition!(state, x -> shu_osher_ic(eq, x))
        result = ssp_rk3!(state, eq, method, tf; cfl=cfl)
        c = run_shu_osher(;
            p=p,
            n_elements=Ne,
            t_final=tf,
            cfl=cfl,
            method=method,
            method_name=opts["method"],
        )
        println("case=shu_osher p=$p ne=$Ne t_final=$tf method=$(opts["method"])")
        println("status=$(c["pass"] ? "ok" : "fail") η=$(c["overshoot"]) δ=$(c["shock_thickness"]) Dex=$(c["excess_dissipation"])")
        println("positivity=$(c["positivity_ok"])")
        if result.status == :ok
            _maybe_write_vtu(opts, state, eq)
        end
        return c["pass"] ? 0 : 1
    elseif case == "advection2d"
        ax, ay = 1.0, 0.5
        eq = LinearAdvection2D(ax, ay)
        ops = build_operators(p)
        mesh = Mesh2D(0.0, 1.0, 0.0, 1.0, Ne, Ne)
        state = allocate_state(mesh, ops, Val(1))
        set_initial_condition!(state, (x, y) -> sin(2π * x) * sin(2π * y))
        result = ssp_rk3!(state, eq, method, t_final; cfl=cfl)
        err = l2_error(
            state,
            (x, y) -> sin(2π * (x - ax * t_final)) * sin(2π * (y - ay * t_final)),
            1,
        )
        println("case=advection2d p=$p ne=$(Ne)x$(Ne) t_final=$t_final L2=$err status=$(result.status)")
        _maybe_write_vtu(opts, state, eq)
        return result.status == :ok ? 0 : 1
    elseif case == "euler2d_wave"
        eq = Euler2D(1.4)
        ops = build_operators(p)
        mesh = Mesh2D(0.0, 1.0, 0.0, 1.0, Ne, Ne)
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
        result = ssp_rk3!(state, eq, method, t_final; cfl=min(cfl, 0.15))
        println(
            "case=euler2d_wave p=$p ne=$(Ne)x$(Ne) t_final=$t_final status=$(result.status) pos=$(positivity_ok(eq, state))",
        )
        _maybe_write_vtu(opts, state, eq)
        return result.status == :ok ? 0 : 1
    elseif case == "euler2d_jump"
        eq = Euler2D(1.4)
        ops = build_operators(p)
        mesh = Mesh2D(
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
        state = allocate_state(mesh, ops, Val(4))
        set_initial_condition!(
            state,
            (x, y) ->
                x < 0.5 ? primitives_to_conserved(eq, 1.0, 0.0, 0.0, 1.0) :
                primitives_to_conserved(eq, 0.125, 0.0, 0.0, 0.1),
        )
        tf = opts["t_final"] == 1.0 ? 0.05 : opts["t_final"]
        result = ssp_rk3!(state, eq, method, tf; cfl=min(cfl, 0.1))
        println(
            "case=euler2d_jump p=$p ne=$(Ne)x$(Ne) t_final=$tf status=$(result.status) pos=$(positivity_ok(eq, state))",
        )
        _maybe_write_vtu(opts, state, eq)
        return result.status == :ok ? 0 : 1
    else
        println(stderr, "Unknown case: $case")
        return 2
    end
end

function cli_invent(args)
    opts = _parse_invent_args(args)
    method = opts["method"]
    baseline = opts["baseline"]
    haskey(METHOD_REGISTRY, method) || error("Unknown method \"$method\". $(describe_methods())")
    haskey(METHOD_REGISTRY, baseline) || error("Unknown baseline \"$baseline\".")
    met, bas, cmp = invent_method(
        method;
        baseline=baseline,
        report_dir=opts["report_dir"],
        δ=opts["delta"],
        vtk_produced=opts["vtk_produced"],
    )
    status = cmp["candidate_status"]
    # Exit 0 for promising/accepted/pass_gates; 1 for rejected
    return status == "rejected" ? 1 : 0
end

function cli_score(args)
    opts = _parse_score_args(args)
    out = isempty(opts["output"]) ? nothing : opts["output"]
    cmp = score_reports(
        opts["method_report"],
        opts["baseline_report"];
        δ=opts["delta"],
        vtk_produced=opts["vtk_produced"],
        out_path=out,
    )
    return cmp["candidate_status"] == "rejected" ? 1 : 0
end

"""
    main_cli(args=ARGS) -> Int

Top-level CLI dispatcher. Returns a process exit code.
"""
function main_cli(args=ARGS)
    if isempty(args) || args[1] in ("-h", "--help", "help")
        _print_usage()
        return isempty(args) ? 2 : 0
    end

    cmd = args[1]
    rest = args[2:end]

    try
        if cmd == "test"
            return cli_test(_parse_test_args(rest))
        elseif cmd == "run"
            return cli_run(rest)
        elseif cmd == "invent"
            return cli_invent(rest)
        elseif cmd == "score"
            return cli_score(rest)
        else
            println(stderr, "Unknown command: $cmd")
            _print_usage()
            return 2
        end
    catch e
        if e isa ArgParseError
            println(stderr, "Argument error: ", e.text)
            return 2
        end
        rethrow()
    end
end
