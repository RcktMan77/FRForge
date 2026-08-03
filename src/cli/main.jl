# CLI entry points for the `frforge` binary.

function _print_usage(io=stderr)
    println(io, "FRForge — high-order Flux Reconstruction laboratory")
    println(io, "Usage: frforge {test|run|invent|score} [options]")
    println(io)
    println(io, "Commands:")
    println(io, "  test    Run verification suite and emit JSON report")
    println(io, "  run     Run a single case")
    println(io, "  invent  Propose/score a capturing method vs baseline (Milestone 6)")
    println(io, "  score   Score two existing reports (Milestone 6)")
    println(io)
    println(io, "Examples:")
    println(io, "  frforge test --suite smoke --report results/smoke/report.json")
    println(io, "  frforge test --suite advection --report results/m1/report.json")
    println(io, "  frforge run --case advection_sine --p 3 --ne 16")
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
        help = "Suite name: smoke | advection | m1 | full"
        default = "advection"
        "--method"
        help = "Capturing method name (null until M4)"
        default = "null"
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
        help = "Case name (advection_sine)"
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
    end
    return parse_args(args, s)
end

"""
    cli_test(opts) -> Int

Run verification suite and write schema v1 JSON report.
Suites: smoke (empty skeleton), advection/m1/full (M1 order + conservation).
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
    elseif suite in ("advection", "m1", "full")
        cases, overall_pass, hard_fails = run_m1_advection_suite()
        diverged = any(c -> get(c, "diverged", false) === true, cases)
        nan_detected = any(c -> get(c, "nan_detected", false) === true, cases)
    else
        println(stderr, "Unknown suite: $suite (use smoke|advection|m1|full)")
        return 2
    end

    report = write_report_skeleton(
        report_path;
        command="test",
        suite=suite,
        method_name=method,
        method_params=Dict{String,Any}(),
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
    println(
        "schema_version=$(report["schema_version"]) overall_pass=$(report["overall_pass"]) n_cases=$(report["summary"]["n_cases"])",
    )
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
        elseif haskey(c, "conservation_residual")
            extra = " cons_res=$(c["conservation_residual"])"
        end
        println("  [$status] $(c["name"])$extra")
    end

    return overall_pass ? 0 : 1
end

function cli_run(args)
    opts = _parse_run_args(args)
    case = opts["case"]
    if case != "advection_sine"
        println(stderr, "Unknown case: $case (M1 supports advection_sine)")
        return 2
    end
    p = opts["p"]
    Ne = opts["ne"]
    a = opts["a"]
    t_final = opts["t_final"]
    cfl = opts["cfl"]

    eq = LinearAdvection1D(a)
    ops = build_operators(p)
    mesh = Mesh1D(0.0, 1.0, Ne; left_bc=PeriodicBC(), right_bc=PeriodicBC())
    state = allocate_state(mesh, ops, Val(1))
    set_initial_condition!(state, x -> sin(2π * x))
    M0 = discrete_mass(state, 1)
    result = ssp_rk3!(state, eq, NullCapturing(), t_final; cfl=cfl)
    MT = discrete_mass(state, 1)
    err = l2_error(state, x -> sin(2π * (x - a * t_final)), 1)

    println("case=advection_sine p=$p ne=$Ne a=$a t_final=$t_final cfl=$cfl")
    println("status=$(result.status) n_steps=$(result.n_steps) t=$(result.t)")
    println("L2_error=$err mass_change=$(MT - M0)")
    return result.status == :ok ? 0 : 1
end

function cli_invent(_opts)
    println(stderr, "frforge invent is not available until Milestone 6.")
    return 2
end

function cli_score(_opts)
    println(stderr, "frforge score is not available until Milestone 6.")
    return 2
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
