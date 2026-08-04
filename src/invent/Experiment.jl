# Invention experiment runner + shared report/JSON artifact helpers.
#
# Layout (invent/):
#   Experiment.jl     — run_method_report, write_report / write_json_pretty, paths
#   Candidate.jl      — classify_candidate vs baseline
#   Invent.jl         — invent_method / score_reports orchestration
#   ExperimentLog.jl  — append/read research/experiment_log.md
#   LogAnalytics.jl   — summary / frontier / lessons views
#   Confirm.jl        — fine-mesh confirm peer of invent
#   Robustness.jl     — scheme-axis matrix
#   Snapshot.jl       — freeze / verify / tables

const DEFAULT_SCORE_MARGIN = 0.02

"""
    run_method_report(method_name; suite=:quant, scheme=DEFAULT_SCHEME) -> Dict

Run the quantitative suite for `method_name` and return a schema v1 report dict
with absolute scores filled.

`suite`: `:quant` / `:m5` (full) or `:light` (CI / robustness reduced).
Invent history must use `DEFAULT_SCHEME` unless a logged re-baseline is recorded.
"""
function run_method_report(
    method_name::AbstractString;
    suite::Symbol = :quant,
    seed = nothing,
    scheme::SchemeConfig = DEFAULT_SCHEME,
)
    # Invent/score path always serial residual (bit-deterministic composite history)
    return with_serial_residual() do
        t0 = time()
        light = suite === :light || suite === :robustness_light
        if suite === :quant || suite === :m5 || light
            cases, overall, hard_fails =
                run_m5_quant_suite(; method_name = method_name, scheme = scheme, light = light)
        else
            error("unsupported suite=$suite (use :quant, :m5, or :light)")
        end
        diverged = any(c -> get(c, "diverged", false) === true, cases)
        nan_detected = any(c -> get(c, "nan_detected", false) === true, cases)
        m = get_capturing_method(method_name)
        report = report_skeleton(;
            command = "invent",
            suite = light ? "light" : String(suite),
            method_name = String(method_name),
            method_params = method_params(m),
            baseline_name = nothing,
            overall_pass = overall,
            diverged = diverged,
            nan_detected = nan_detected,
            wall_time_sec = time() - t0,
            hard_gate_failures = hard_fails,
            cases = cases,
            fill_scores = true,
            scheme = scheme,
        )
        if seed !== nothing
            report["method_params"]["seed"] = seed
        end
        return report
    end
end

# ---------------------------------------------------------------------------
# Shared report / JSON artifacts (used by invent, confirm, robustness, CLI)
# ---------------------------------------------------------------------------

"""
    write_json_pretty(path, obj) -> path

Write any JSON-serializable value pretty-printed with a trailing newline.
Creates parent directories. Does **not** validate report schema keys.
"""
function write_json_pretty(path::AbstractString, obj)
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        JSON.print(io, obj, 2)
        println(io)
    end
    return path
end

"""
    write_report(path, report) -> report

Validate schema keys, then write pretty JSON (via `write_json_pretty`).
"""
function write_report(path::AbstractString, report::AbstractDict)
    errs = validate_report_keys(report)
    isempty(errs) || error("report validation failed: $(join(errs, "; "))")
    write_json_pretty(path, report)
    return report
end

"""
    report_trio_paths(report_dir, method, baseline; tag="") -> NamedTuple

Standard method / baseline / compare JSON paths.

- invent: `tag=""` → `method_<m>.json`, `baseline_<b>.json`, `compare_<m>_vs_<b>.json`
- confirm: `tag=preset` → same with `_<preset>` suffix before `.json`
"""
function report_trio_paths(
    report_dir::AbstractString,
    method::AbstractString,
    baseline::AbstractString;
    tag::AbstractString = "",
)
    suffix = isempty(tag) ? "" : "_$(tag)"
    return (
        method = joinpath(report_dir, "method_$(method)$(suffix).json"),
        baseline = joinpath(report_dir, "baseline_$(baseline)$(suffix).json"),
        compare = joinpath(report_dir, "compare_$(method)_vs_$(baseline)$(suffix).json"),
    )
end

"""
Stamp invent/confirm bookkeeping fields on a report dict (in-place).
"""
function stamp_workflow_report!(
    report::AbstractDict;
    command::AbstractString,
    baseline_name = nothing,
    scheme::SchemeConfig = DEFAULT_SCHEME,
    extra::AbstractDict = Dict{String,Any}(),
)
    report["command"] = command
    report["baseline_name"] = baseline_name
    report["scheme"] = scheme_dict(scheme)
    for (k, v) in extra
        report[k] = v
    end
    return report
end

"""Print the usual three report paths to stdout."""
function print_report_trio(met_path, bas_path, cmp_path)
    println("Reports: $met_path")
    println("         $bas_path")
    println("         $cmp_path")
    return nothing
end

"""
    resolve_log_paths(log_path, yaml_path) -> (log_path, yaml_path_or_nothing)

Default experiment-log paths; empty-string `yaml_path` disables the YAML index.
"""
function resolve_log_paths(log_path, yaml_path)
    lp = something(log_path, default_experiment_log_path())
    yp = yaml_path === nothing ? default_experiment_log_yaml_path() : yaml_path
    yp_use = (yp isa AbstractString && isempty(yp)) ? nothing : yp
    return lp, yp_use
end

"""
    report_artifact_dict(met_path, bas_path, cmp_path; extra...) -> Dict

Standard invent/confirm artifact map for experiment-log entries.
"""
function report_artifact_dict(met_path, bas_path, cmp_path; extra::AbstractDict=Dict{String,Any}())
    arts = Dict{String,Any}(
        "method_report" => met_path,
        "baseline_report" => bas_path,
        "compare" => cmp_path,
    )
    for (k, v) in pairs(extra)
        arts[String(k)] = v
    end
    return arts
end

"""
    maybe_append_workflow_log!(; append_log, entry_builder, log_path, yaml_path) -> entry|nothing

Call `entry_builder(log_path, yaml_path)` when `append_log` is true, print the
appended id, and warn on incomplete narrative when present.
`entry_builder` should return the entry dict (with optional `narrative_complete`).
"""
function maybe_append_workflow_log!(;
    append_log::Bool,
    entry_builder::Function,
    log_path=nothing,
    yaml_path=nothing,
)
    append_log || return nothing
    lp, yp = resolve_log_paths(log_path, yaml_path)
    entry = entry_builder(lp, yp)
    println("Experiment log appended: $(entry["id"])  →  $lp")
    if get(entry, "narrative_complete", true) === false
        println(
            "  WARNING: hypothesis/lessons are placeholders — required for promising+ before shortlist.",
        )
    end
    return entry
end
