# Top-level invent / score orchestration.

"""
    invent_method(method_name; baseline="persson_av", report_dir="results/invent", δ=0.02, vtk_produced=false)

Run quant suite for method and baseline, write JSON reports, classify candidate.
Returns (method_report, baseline_report, comparison_dict).
"""
function invent_method(
    method_name::AbstractString;
    baseline::AbstractString="persson_av",
    report_dir::AbstractString="results/invent",
    δ::Real=DEFAULT_SCORE_MARGIN,
    vtk_produced::Bool=false,
    seed=nothing,
)
    mkpath(report_dir)
    bas_path = joinpath(report_dir, "baseline_$(baseline).json")
    met_path = joinpath(report_dir, "method_$(method_name).json")
    cmp_path = joinpath(report_dir, "compare_$(method_name)_vs_$(baseline).json")

    println("Running baseline suite: $baseline ...")
    bas = run_method_report(baseline; seed=seed)
    bas["command"] = "invent"
    bas["baseline_name"] = nothing
    write_report(bas_path, bas)

    println("Running method suite: $method_name ...")
    met = run_method_report(method_name; seed=seed)
    met["command"] = "invent"
    met["baseline_name"] = baseline
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
