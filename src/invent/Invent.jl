# Top-level invent / score orchestration (coarse quant + classify).
# Fine-mesh path lives in Confirm.jl; scheme matrix in Robustness.jl.

"""Attach classify_candidate fields onto the method report dict (in-place)."""
function _attach_candidate_fields!(met::AbstractDict, cmp::AbstractDict)
    met["candidate_status"] = cmp["candidate_status"]
    met["baseline_composite"] = cmp["baseline_composite"]
    met["composite_margin"] = cmp["composite_margin"]
    met["score_margin_threshold"] = cmp["score_margin_threshold"]
    met["vtk_produced"] = cmp["vtk_produced"]
    met["tradeoff_ok"] = cmp["tradeoff_ok"]
    met["summary"]["tradeoff_notes"] = cmp["tradeoff_notes"]
    met["summary"]["scores"] = cmp["absolute_scores"]
    return met
end

function _invent_append_if_requested!(
    method_name,
    met,
    bas,
    cmp;
    append_log::Bool,
    log_path,
    yaml_path,
    met_path,
    bas_path,
    cmp_path,
    hypothesis,
    lessons,
    strengths,
    weaknesses,
    git_ref,
)
    return maybe_append_workflow_log!(;
        append_log=append_log,
        log_path=log_path,
        yaml_path=yaml_path,
        entry_builder=(lp, yp) -> invent_append_log!(
            method_name,
            met,
            bas,
            cmp;
            log_path=lp,
            yaml_path=yp,
            artifacts=report_artifact_dict(met_path, bas_path, cmp_path),
            hypothesis=hypothesis,
            lessons=lessons,
            strengths=strengths,
            weaknesses=weaknesses,
            git_ref=git_ref,
        ),
    )
end

"""
    invent_method(method_name; baseline="persson_av", report_dir="results/invent",
                  δ=0.02, vtk_produced=false, append_log=true, ...)

Run quant suite for method and baseline, write JSON reports, classify candidate.
By default appends a stub entry to `research/experiment_log.md` (authoritative memory).

**Scheme:** invent always uses the frozen default scheme (GL + Rusanov + SSP-RK3)
for composite-score history unless a logged re-baseline is recorded.

Returns (method_report, baseline_report, comparison_dict).
"""
function invent_method(
    method_name::AbstractString;
    baseline::AbstractString = "persson_av",
    report_dir::AbstractString = "results/invent",
    δ::Real = DEFAULT_SCORE_MARGIN,
    vtk_produced::Bool = false,
    seed = nothing,
    append_log::Bool = true,
    log_path::Union{Nothing,AbstractString} = nothing,
    yaml_path::Union{Nothing,AbstractString} = nothing,
    hypothesis::AbstractString = "",
    lessons::AbstractString = "",
    strengths::AbstractString = "",
    weaknesses::AbstractString = "",
    git_ref::AbstractString = "",
)
    method_name = require_registered_method(method_name)
    baseline = require_registered_method(baseline; role = "baseline")

    mkpath(report_dir)
    paths = report_trio_paths(report_dir, method_name, baseline)
    bas_path, met_path, cmp_path = paths.baseline, paths.method, paths.compare

    println("Reminder: read research/experiment_log.md before proposing new methods.")
    println(
        "Frozen invent scheme: ",
        FROZEN_INVENT_SCHEME.points,
        " + ",
        FROZEN_INVENT_SCHEME.flux,
        " + ",
        FROZEN_INVENT_SCHEME.time,
    )

    println("Running baseline suite: $baseline ...")
    bas = run_method_report(baseline; seed = seed)
    stamp_workflow_report!(
        bas;
        command = "invent",
        baseline_name = nothing,
        scheme = DEFAULT_SCHEME,
    )
    write_report(bas_path, bas)

    println("Running method suite: $method_name ...")
    met = run_method_report(method_name; seed = seed)
    stamp_workflow_report!(
        met;
        command = "invent",
        baseline_name = baseline,
        scheme = DEFAULT_SCHEME,
    )
    write_report(met_path, met)

    cmp = classify_candidate(met, bas; δ = δ, vtk_produced = vtk_produced)
    _attach_candidate_fields!(met, cmp)
    write_report(met_path, met)
    write_json_pretty(cmp_path, cmp)

    print_candidate_summary(method_name, cmp)
    print_report_trio(met_path, bas_path, cmp_path)

    _invent_append_if_requested!(
        method_name,
        met,
        bas,
        cmp;
        append_log = append_log,
        log_path = log_path,
        yaml_path = yaml_path,
        met_path = met_path,
        bas_path = bas_path,
        cmp_path = cmp_path,
        hypothesis = hypothesis,
        lessons = lessons,
        strengths = strengths,
        weaknesses = weaknesses,
        git_ref = git_ref,
    )

    return met, bas, cmp
end

"""
    score_reports(method_path, baseline_path; δ=0.02, vtk_produced=false) -> Dict

Load two existing JSON reports and classify the method vs baseline.
"""
function score_reports(
    method_path::AbstractString,
    baseline_path::AbstractString;
    δ::Real = DEFAULT_SCORE_MARGIN,
    vtk_produced::Bool = false,
    out_path::Union{Nothing,AbstractString} = nothing,
)
    met = load_report(method_path)
    bas = load_report(baseline_path)
    cmp = classify_candidate(met, bas; δ = δ, vtk_produced = vtk_produced)
    print_candidate_summary(get(met, "method_name", "method"), cmp)
    if out_path !== nothing
        write_json_pretty(out_path, cmp)
        println("Wrote $out_path")
    end
    return cmp
end
