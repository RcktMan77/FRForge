# Candidate status classification vs baseline (invent / score).

"""
    score_suite_relative(method_cases, baseline_cases) -> Dict

Relative S_i maps (formula version 1) when a baseline report is available.
"""
function score_suite_relative(method_cases::AbstractVector, baseline_cases::AbstractVector)
    # Order: still binary on method's smooth_order cases
    order_cases = filter(c -> get(c, "case_type", "") == "smooth_order", method_cases)
    S_order = if isempty(order_cases)
        1.0
    else
        all(c -> get(c, "order_pass", false) === true, order_cases) ? 1.0 : 0.0
    end

    m_disc = filter(c -> get(c, "case_type", "") == "discontinuous", method_cases)
    b_disc = filter(c -> get(c, "case_type", "") == "discontinuous", baseline_cases)

    # Pair by base case name prefix (sod_ / shu_osher_) when possible
    function match_baseline(mc)
        mname = get(mc, "name", "")
        for bc in b_disc
            bname = get(bc, "name", "")
            # strip method suffix after last _
            if startswith(mname, "sod") && startswith(bname, "sod")
                return bc
            end
            if startswith(mname, "shu_osher") && startswith(bname, "shu_osher")
                return bc
            end
        end
        return isempty(b_disc) ? nothing : b_disc[1]
    end

    diss_scores = Float64[]
    shock_scores = Float64[]
    for mc in m_disc
        bc = match_baseline(mc)
        ηM = float(get(mc, "overshoot", 0.0))
        DexM = get(mc, "excess_dissipation", nothing)
        δM = get(mc, "shock_thickness", nothing)
        if bc === nothing
            # fall back absolute
            if DexM === nothing
                push!(diss_scores, clip01(1 - ηM / ETA_REF_ABS))
            else
                push!(diss_scores, clip01(1 - float(DexM) / D_REF_ABS))
            end
            sδ =
                δM === nothing || !isfinite(float(δM)) ? 1.0 :
                clip01(1 - float(δM) / DELTA_REF_ABS)
            push!(shock_scores, 0.5 * sδ + 0.5 * clip01(1 - ηM / ETA_REF_ABS))
            continue
        end
        ηB = float(get(bc, "overshoot", 0.0))
        DexB = get(bc, "excess_dissipation", nothing)
        δB = get(bc, "shock_thickness", nothing)
        if DexM === nothing || DexB === nothing
            push!(diss_scores, clip01(ηB / (ηM + SCORE_EPS)))
        else
            push!(diss_scores, clip01(float(DexB) / (float(DexM) + SCORE_EPS)))
        end
        sδ =
            if δM === nothing || δB === nothing || !isfinite(float(δM)) || !isfinite(float(δB))
                1.0
            else
                clip01(float(δB) / (float(δM) + SCORE_EPS))
            end
        sη = clip01(ηB / (ηM + SCORE_EPS))
        push!(shock_scores, 0.5 * sδ + 0.5 * sη)
    end
    S_diss = isempty(diss_scores) ? 1.0 : sum(diss_scores) / length(diss_scores)
    S_shock = isempty(shock_scores) ? 1.0 : sum(shock_scores) / length(shock_scores)

    robust_fail = any(
        c ->
            get(c, "diverged", false) ||
            get(c, "nan_detected", false) ||
            get(c, "positivity_ok", true) === false,
        method_cases,
    )
    S_robust = robust_fail ? 0.0 : 1.0

    w = DEFAULT_SCORING_WEIGHTS
    composite =
        w["order_preservation"] * S_order +
        w["dissipation"] * S_diss +
        w["shock_quality"] * S_shock +
        w["robustness"] * S_robust

    return Dict{String,Float64}(
        "order_preservation" => S_order,
        "dissipation" => S_diss,
        "shock_quality" => S_shock,
        "robustness" => S_robust,
        "composite" => composite,
    )
end

"""
    tradeoff_ok(method_scores, baseline_scores) -> (Bool, String)

Order-vs-dissipation trade-off: method should not lose order, and should improve
dissipation and/or shock quality (or win composite without order loss).
"""
function tradeoff_ok(method_scores::AbstractDict, baseline_scores::AbstractDict)
    So_m = float(method_scores["order_preservation"])
    So_b = float(baseline_scores["order_preservation"])
    Sd_m = float(method_scores["dissipation"])
    Sd_b = float(baseline_scores["dissipation"])
    Ss_m = float(method_scores["shock_quality"])
    Ss_b = float(baseline_scores["shock_quality"])
    Sc_m = float(method_scores["composite"])
    Sc_b = float(baseline_scores["composite"])

    order_ok = So_m >= So_b - 1e-12
    improved_diss = Sd_m > Sd_b + 1e-12
    improved_shock = Ss_m > Ss_b + 1e-12
    composite_win = Sc_m > Sc_b + 1e-12

    ok = order_ok && (improved_diss || improved_shock || composite_win)
    notes = String[]
    push!(notes, order_ok ? "order_preserved" : "order_regressed")
    improved_diss && push!(notes, "dissipation_improved")
    improved_shock && push!(notes, "shock_improved")
    composite_win && push!(notes, "composite_improved")
    return ok, join(notes, ",")
end

"""
    classify_candidate(method_report, baseline_report; δ=0.02, vtk_produced=false) -> Dict

Return candidate_status and comparison fields per design.
"""
function classify_candidate(
    method_report::AbstractDict,
    baseline_report::AbstractDict;
    δ::Real = DEFAULT_SCORE_MARGIN,
    vtk_produced::Bool = false,
)
    overall = get(method_report, "overall_pass", false) === true
    m_scores = method_report["summary"]["scores"]
    b_scores = baseline_report["summary"]["scores"]
    # Prefer absolute scores already on reports for composite comparison;
    # also compute relative scores for invent ranking display.
    rel = score_suite_relative(method_report["cases"], baseline_report["cases"])
    # Use absolute composites for margin (stable, comparable)
    c_m = float(m_scores["composite"])
    c_b = float(b_scores["composite"])
    margin = c_m - c_b
    t_ok, t_notes = tradeoff_ok(m_scores, b_scores)

    status =
        if !overall || get(method_report, "diverged", false) ||
           get(method_report, "nan_detected", false)
            "rejected"
        elseif margin >= float(δ) && t_ok
            vtk_produced ? "accepted_candidate" : "promising"
        else
            "pass_gates"
        end

    return Dict{String,Any}(
        "candidate_status" => status,
        "baseline_name" => get(baseline_report, "method_name", nothing),
        "baseline_composite" => c_b,
        "composite" => c_m,
        "composite_margin" => margin,
        "score_margin_threshold" => float(δ),
        "vtk_produced" => vtk_produced,
        "tradeoff_ok" => t_ok,
        "tradeoff_notes" => t_notes,
        "relative_scores" => rel,
        "absolute_scores" => m_scores,
        "baseline_absolute_scores" => b_scores,
    )
end

"""
    print_candidate_summary(method_name, cmp; io=stdout)

Human-readable invent/score summary required by design.
"""
function print_candidate_summary(method_name::AbstractString, cmp::AbstractDict; io::IO = stdout)
    status = cmp["candidate_status"]
    baseline = something(cmp["baseline_name"], "none")
    println(io, "Method: $method_name    status: $status")
    println(
        io,
        "overall_pass: see report   composite: $(round(cmp["composite"], digits=4))  baseline($baseline): $(round(cmp["baseline_composite"], digits=4))  margin: $(cmp["composite_margin"] >= 0 ? "+" : "")$(round(cmp["composite_margin"], digits=4))",
    )
    println(
        io,
        "tradeoff_ok: $(cmp["tradeoff_ok"])    vtk_produced: $(cmp["vtk_produced"])    tradeoff: $(cmp["tradeoff_notes"])",
    )
    println(io, "δ_score threshold: $(cmp["score_margin_threshold"])")
    return nothing
end
