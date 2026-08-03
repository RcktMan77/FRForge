# CLI entry points for the `frforge` binary.

function _print_usage(io=stderr)
    println(io, "FRForge — high-order Flux Reconstruction laboratory")
    println(io, "Usage: frforge {test|run|invent|score} [options]")
    println(io)
    println(io, "Commands:")
    println(io, "  test    Run verification suite (or M0 smoke) and emit JSON report")
    println(io, "  run     Run a single case (available from Milestone 1)")
    println(io, "  invent  Propose/score a capturing method vs baseline (Milestone 6)")
    println(io, "  score   Score two existing reports (Milestone 6)")
    println(io)
    println(io, "Examples:")
    println(io, "  frforge test --report results/smoke/report.json")
    println(io, "  frforge test --suite smoke --report out.json")
end

function _parse_test_args(args)
    s = ArgParseSettings(
        description = "Run FRForge verification and write a JSON report.",
        prog = "frforge test",
    )
    @add_arg_table! s begin
        "--report", "-r"
            help = "Path for the JSON report"
            default = "results/report.json"
            dest_name = "report"
        "--suite"
            help = "Suite name (smoke | full | ...)"
            default = "smoke"
        "--method"
            help = "Capturing method name (null until M4)"
            default = "null"
    end
    return parse_args(args, s)
end

"""
    cli_test(opts) -> Int

Milestone 0: write a valid JSON report skeleton and exit 0 if valid.
Later milestones append real cases under the same command.
"""
function cli_test(opts::AbstractDict)
    t0 = time()
    report_path = opts["report"]
    suite = opts["suite"]
    method = opts["method"]

    # M0: smoke suite produces an empty-but-valid skeleton.
    # Later milestones will populate `cases` with real verification results.
    report = write_report_skeleton(
        report_path;
        command = "test",
        suite = suite,
        method_name = method,
        method_params = Dict{String,Any}(),
        baseline_name = nothing,
        overall_pass = true,
        diverged = false,
        nan_detected = false,
        wall_time_sec = time() - t0,
        hard_gate_failures = String[],
        cases = Any[],
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
    println("schema_version=$(report["schema_version"]) overall_pass=$(report["overall_pass"])")
    return 0
end

function cli_run(_opts)
    println(stderr, "frforge run is not available until Milestone 1.")
    return 2
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
