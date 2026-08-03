# JSON verification report writer and key-list validator.

"""
    package_version() -> String

Return the package version string from `Project.toml` (best-effort).
"""
function package_version()
    try
        proj = joinpath(@__DIR__, "..", "..", "Project.toml")
        for line in eachline(proj)
            m = match(r"^version\s*=\s*\"([^\"]+)\"", line)
            m !== nothing && return String(m.captures[1])
        end
    catch
    end
    return "0.1.0"
end

"""
    git_commit_short() -> String

Best-effort short git SHA of the working tree, or `"unknown"`.
"""
function git_commit_short()
    try
        root = joinpath(@__DIR__, "..", "..")
        sha = read(Cmd(`git -C $root rev-parse --short HEAD`; ignorestatus=true), String)
        sha = strip(sha)
        return isempty(sha) ? "unknown" : sha
    catch
        return "unknown"
    end
end

"""
    empty_scores() -> Dict

Placeholder score components (null until M5/M6 fill them).
"""
function empty_scores()
    return Dict{String,Any}(
        "order_preservation" => nothing,
        "dissipation" => nothing,
        "shock_quality" => nothing,
        "robustness" => nothing,
        "composite" => nothing,
    )
end

"""
    report_skeleton(; kwargs...) -> Dict{String,Any}

Build a schema v1 report skeleton with no cases.
Valid for M0; later milestones append cases and fill scores.
"""
function report_skeleton(;
    command::AbstractString = "test",
    suite::AbstractString = "smoke",
    method_name::AbstractString = "null",
    method_params = Dict{String,Any}(),
    baseline_name = nothing,
    overall_pass::Bool = true,
    diverged::Bool = false,
    nan_detected::Bool = false,
    wall_time_sec::Real = 0.0,
    hard_gate_failures = String[],
    cases = Any[],
    fill_scores::Bool = true,
)
    n_cases = length(cases)
    n_passed = count(c -> get(c, "pass", false) === true, cases)
    n_failed = n_cases - n_passed

    scores = empty_scores()
    if fill_scores && n_cases > 0
        scores = score_suite_absolute(cases)
    end

    return Dict{String,Any}(
        "schema_version" => SCHEMA_VERSION,
        "package" => "FRForge",
        "package_version" => package_version(),
        "git_commit" => git_commit_short(),
        "timestamp_utc" => Dates.format(Dates.now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS") * "Z",
        "julia_version" => string(VERSION),
        "command" => String(command),
        "suite" => String(suite),
        "method_name" => String(method_name),
        "method_params" => method_params,
        "baseline_name" => baseline_name,
        "overall_pass" => overall_pass,
        "diverged" => diverged,
        "nan_detected" => nan_detected,
        "wall_time_sec" => Float64(wall_time_sec),
        "scoring_weights" => deepcopy(DEFAULT_SCORING_WEIGHTS),
        "scoring_formula_version" => SCORING_FORMULA_VERSION,
        "hard_gate_failures" => collect(hard_gate_failures),
        "cases" => collect(cases),
        "summary" => Dict{String,Any}(
            "n_cases" => n_cases,
            "n_passed" => n_passed,
            "n_failed" => n_failed,
            "scores" => scores,
        ),
    )
end

"""
    validate_report_keys(report) -> Vector{String}

Return a list of missing/invalid key messages. Empty means the report
satisfies the schema v1 key contract (types lightly checked).
"""
function validate_report_keys(report)::Vector{String}
    errors = String[]
    if !(report isa AbstractDict)
        push!(errors, "report is not a dictionary")
        return errors
    end

    for k in REQUIRED_TOP_LEVEL_KEYS
        if !haskey(report, k)
            push!(errors, "missing top-level key: $k")
        end
    end
    isempty(errors) || return errors

    if report["schema_version"] != SCHEMA_VERSION
        push!(errors, "schema_version must be $SCHEMA_VERSION, got $(report["schema_version"])")
    end
    if report["scoring_formula_version"] != SCORING_FORMULA_VERSION
        push!(
            errors,
            "scoring_formula_version must be $SCORING_FORMULA_VERSION, got $(report["scoring_formula_version"])",
        )
    end
    if haskey(report, "overall_score")
        push!(errors, "overall_score is forbidden; use summary.scores.composite only")
    end

    summary = report["summary"]
    if !(summary isa AbstractDict)
        push!(errors, "summary must be an object")
        return errors
    end
    for k in REQUIRED_SUMMARY_KEYS
        if !haskey(summary, k)
            push!(errors, "missing summary key: $k")
        end
    end

    scores = get(summary, "scores", nothing)
    if scores isa AbstractDict
        for k in REQUIRED_SCORE_KEYS
            if !haskey(scores, k)
                push!(errors, "missing summary.scores key: $k")
            end
        end
    else
        push!(errors, "summary.scores must be an object")
    end

    cases = report["cases"]
    if !(cases isa AbstractVector)
        push!(errors, "cases must be an array")
        return errors
    end
    for (i, case) in enumerate(cases)
        if !(case isa AbstractDict)
            push!(errors, "cases[$i] must be an object")
            continue
        end
        for k in REQUIRED_CASE_KEYS
            if !haskey(case, k)
                push!(errors, "cases[$i] missing key: $k")
            end
        end
    end

    return errors
end

"""
    write_report_skeleton(path; kwargs...) -> Dict{String,Any}

Create parent directories if needed, write a valid schema v1 skeleton JSON
report, and return the report dict.
"""
function write_report_skeleton(path::AbstractString; kwargs...)
    report = report_skeleton(; kwargs...)
    errs = validate_report_keys(report)
    if !isempty(errs)
        error("internal skeleton failed validation: $(join(errs, "; "))")
    end
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        JSON.print(io, report, 2)
        println(io)
    end
    return report
end

"""
    load_report(path) -> Dict

Parse a JSON report file into a Julia dictionary.
"""
function load_report(path::AbstractString)
    return JSON.parsefile(path)
end
