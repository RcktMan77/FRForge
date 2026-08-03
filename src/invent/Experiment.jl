# Invention experiment runner: run quant suite for a method and optional baseline.

const DEFAULT_SCORE_MARGIN = 0.02

"""
    run_method_report(method_name; suite=:quant, report_path=nothing) -> Dict

Run the quantitative suite for `method_name` and return a schema v1 report dict
with absolute scores filled.
"""
function run_method_report(
    method_name::AbstractString;
    suite::Symbol=:quant,
    seed=nothing,
)
    t0 = time()
    if suite === :quant || suite === :m5
        cases, overall, hard_fails = run_m5_quant_suite(; method_name=method_name)
    else
        error("invent currently supports suite=quant only (got $suite)")
    end
    diverged = any(c -> get(c, "diverged", false) === true, cases)
    nan_detected = any(c -> get(c, "nan_detected", false) === true, cases)
    m = get_capturing_method(method_name)
    report = report_skeleton(;
        command="invent",
        suite=String(suite),
        method_name=String(method_name),
        method_params=method_params(m),
        baseline_name=nothing,
        overall_pass=overall,
        diverged=diverged,
        nan_detected=nan_detected,
        wall_time_sec=time() - t0,
        hard_gate_failures=hard_fails,
        cases=cases,
        fill_scores=true,
    )
    if seed !== nothing
        report["method_params"]["seed"] = seed
    end
    return report
end

"""
    write_report(path, report) -> report

Ensure parent dirs exist and write pretty JSON.
"""
function write_report(path::AbstractString, report::AbstractDict)
    errs = validate_report_keys(report)
    isempty(errs) || error("report validation failed: $(join(errs, "; "))")
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        JSON.print(io, report, 2)
        println(io)
    end
    return report
end
