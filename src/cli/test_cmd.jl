# `frforge test` — verification suite runner (serial residual).

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
        help = CLI_TEST_SUITE_HELP
        default = "advection"
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
    end
    return parse_args(args, s)
end

"""Build SchemeConfig from CLI option dict (points/flux/time keys)."""
_scheme_from_opts(opts::AbstractDict) = scheme_from_cli_opts(opts)

"""
    cli_test(opts) -> Int

Run verification suite and write schema v1 JSON report.
See `CLI_TEST_SUITE_HELP` for suite names. Always forces serial residual.
"""
function cli_test(opts::AbstractDict)
    # Force serial residual for all CLI test suites (CI / score history determinism)
    return with_serial_residual() do
        return _cli_test_body(opts)
    end
end

"""Case flags derived from a verification case list."""
function _case_flags(cases)
    diverged = any(c -> get(c, "diverged", false) === true, cases)
    nan_detected = any(c -> get(c, "nan_detected", false) === true, cases)
    return diverged, nan_detected
end

"""
Dispatch `frforge test --suite` to suite runners.

Returns `(cases, overall_pass, hard_fails, method_name, ok)` where `ok=false`
means the suite name was unknown (caller should exit 2).
"""
function _dispatch_test_suite(suite::AbstractString, method::AbstractString)
    hard_fails = String[]
    if suite == "smoke"
        return Any[], true, hard_fails, method, true
    elseif suite in ("advection", "m1")
        cases, overall_pass, hard_fails = run_m1_advection_suite()
        return cases, overall_pass, hard_fails, method, true
    elseif suite in ("burgers", "m2")
        cases, overall_pass, hard_fails = run_m2_burgers_suite()
        return cases, overall_pass, hard_fails, method, true
    elseif suite in ("euler", "m3")
        cases, overall_pass, hard_fails = run_m3_euler_suite()
        return cases, overall_pass, hard_fails, method, true
    elseif suite in ("capturing", "m4")
        cases, overall_pass, hard_fails = run_m4_capturing_suite()
        return cases, overall_pass, hard_fails, "persson_av", true
    elseif suite in ("quant", "m5")
        cases, overall_pass, hard_fails = run_m5_quant_suite(; method_name=method)
        return cases, overall_pass, hard_fails, method, true
    elseif suite in ("2d", "m8")
        cases, overall_pass, hard_fails, _, _ = run_m8_2d_suite()
        return cases, overall_pass, hard_fails, method, true
    elseif suite in ("2d_capturing", "p31")
        cases, overall_pass, hard_fails = run_p31_2d_capturing_suite()
        return cases, overall_pass, hard_fails, "persson_av", true
    elseif suite in ("curved", "p32")
        cases, overall_pass, hard_fails = run_p32_curved_suite()
        return cases, overall_pass, hard_fails, "null", true
    elseif suite in ("benchmarks", "p33a", "riemann", "vortex")
        cases, overall_pass, hard_fails = run_p33a_benchmark_suite()
        return cases, overall_pass, hard_fails, "mixed", true
    elseif suite in ("optional2d", "p33b", "dmr", "ffs")
        cases, overall_pass, hard_fails = run_p33b_optional_suite()
        return cases, overall_pass, hard_fails, "persson_av", true
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
        return cases, overall_pass, hard_fails, method, true
    else
        return Any[], false, hard_fails, method, false
    end
end

function _cli_test_body(opts::AbstractDict)
    t0 = time()
    report_path = opts["report"]
    suite = opts["suite"]
    method = opts["method"]
    scheme = _scheme_from_opts(opts)
    # Required CI / invent history use DEFAULT_SCHEME only; non-default is for local exploration
    if scheme != DEFAULT_SCHEME
        println(
            "Note: non-default scheme $(scheme_dict(scheme)) — invent composite history remains frozen at DEFAULT_SCHEME.",
        )
    end

    cases, overall_pass, hard_fails, method, suite_ok =
        _dispatch_test_suite(suite, method)
    if !suite_ok
        println(stderr, "Unknown suite: $suite")
        println(stderr, "  ", CLI_TEST_SUITE_HELP)
        return 2
    end
    diverged, nan_detected = _case_flags(cases)

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
        scheme=scheme,
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
