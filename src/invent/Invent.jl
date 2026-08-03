# Top-level invent / score orchestration.

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
    baseline::AbstractString="persson_av",
    report_dir::AbstractString="results/invent",
    δ::Real=DEFAULT_SCORE_MARGIN,
    vtk_produced::Bool=false,
    seed=nothing,
    append_log::Bool=true,
    log_path::Union{Nothing,AbstractString}=nothing,
    yaml_path::Union{Nothing,AbstractString}=nothing,
    hypothesis::AbstractString="",
    lessons::AbstractString="",
    strengths::AbstractString="",
    weaknesses::AbstractString="",
    git_ref::AbstractString="",
)
    mkpath(report_dir)
    bas_path = joinpath(report_dir, "baseline_$(baseline).json")
    met_path = joinpath(report_dir, "method_$(method_name).json")
    cmp_path = joinpath(report_dir, "compare_$(method_name)_vs_$(baseline).json")

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
    bas = run_method_report(baseline; seed=seed)
    bas["command"] = "invent"
    bas["baseline_name"] = nothing
    bas["scheme"] = scheme_dict(DEFAULT_SCHEME)
    write_report(bas_path, bas)

    println("Running method suite: $method_name ...")
    met = run_method_report(method_name; seed=seed)
    met["command"] = "invent"
    met["baseline_name"] = baseline
    met["scheme"] = scheme_dict(DEFAULT_SCHEME)
    write_report(met_path, met)

    cmp = classify_candidate(met, bas; δ=δ, vtk_produced=vtk_produced)
    # Attach comparison onto method report for agent convenience
    met["candidate_status"] = cmp["candidate_status"]
    met["baseline_composite"] = cmp["baseline_composite"]
    met["composite_margin"] = cmp["composite_margin"]
    met["score_margin_threshold"] = cmp["score_margin_threshold"]
    met["vtk_produced"] = cmp["vtk_produced"]
    met["tradeoff_ok"] = cmp["tradeoff_ok"]
    met["summary"]["tradeoff_notes"] = cmp["tradeoff_notes"]
    met["summary"]["scores"] = cmp["absolute_scores"]
    # re-write method report with invent fields
    write_report(met_path, met)

    open(cmp_path, "w") do io
        JSON.print(io, cmp, 2)
        println(io)
    end

    print_candidate_summary(method_name, cmp)
    println("Reports: $met_path")
    println("         $bas_path")
    println("         $cmp_path")

    if append_log
        lp = something(log_path, default_experiment_log_path())
        yp = yaml_path === nothing ? default_experiment_log_yaml_path() : yaml_path
        # empty string yaml_path disables YAML index update
        yp_use = (yp isa AbstractString && isempty(yp)) ? nothing : yp
        arts = Dict{String,Any}(
            "method_report" => met_path,
            "baseline_report" => bas_path,
            "compare" => cmp_path,
        )
        entry = invent_append_log!(
            method_name,
            met,
            bas,
            cmp;
            log_path=lp,
            yaml_path=yp_use,
            artifacts=arts,
            hypothesis=hypothesis,
            lessons=lessons,
            strengths=strengths,
            weaknesses=weaknesses,
            git_ref=git_ref,
        )
        println("Experiment log appended: $(entry["id"])  →  $lp")
        if get(entry, "narrative_complete", true) === false
            println(
                "  WARNING: hypothesis/lessons are placeholders — required for promising+ before shortlist.",
            )
        end
    end

    return met, bas, cmp
end

"""
    score_reports(method_path, baseline_path; δ=0.02, vtk_produced=false) -> Dict

Load two existing JSON reports and classify the method vs baseline.
"""
function score_reports(
    method_path::AbstractString,
    baseline_path::AbstractString;
    δ::Real=DEFAULT_SCORE_MARGIN,
    vtk_produced::Bool=false,
    out_path::Union{Nothing,AbstractString}=nothing,
)
    met = load_report(method_path)
    bas = load_report(baseline_path)
    cmp = classify_candidate(met, bas; δ=δ, vtk_produced=vtk_produced)
    print_candidate_summary(get(met, "method_name", "method"), cmp)
    if out_path !== nothing
        open(out_path, "w") do io
            JSON.print(io, cmp, 2)
            println(io)
        end
        println("Wrote $out_path")
    end
    return cmp
end
