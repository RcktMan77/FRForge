# Absolute and relative scoring maps (formula version 1). See docs/design.md.

const D_REF_ABS = 0.2
const DELTA_REF_ABS = 8.0   # shock thickness in SP spacings
const ETA_REF_ABS = 0.1
const SCORE_EPS = 1e-30

clip01(x) = clamp(float(x), 0.0, 1.0)

"""
    score_suite_absolute(cases) -> Dict{String,Float64}

Compute S_order, S_diss, S_shock, S_robust, composite from a vector of case dicts
using absolute maps (no baseline). Used by quant/test reports and invent absolute scores.
"""
function score_suite_absolute(cases::AbstractVector)
    # Order
    order_cases = filter(c -> get(c, "case_type", "") == "smooth_order", cases)
    if isempty(order_cases)
        S_order = 1.0
    else
        S_order = all(c -> get(c, "order_pass", false) === true, order_cases) ? 1.0 : 0.0
    end

    # Dissipation & shock from discontinuous cases
    disc = filter(c -> get(c, "case_type", "") == "discontinuous", cases)
    diss_scores = Float64[]
    shock_scores = Float64[]
    for c in disc
        Dex = get(c, "excess_dissipation", nothing)
        η = float(get(c, "overshoot", 0.0))
        if Dex === nothing || Dex === missing
            push!(diss_scores, clip01(1 - η / ETA_REF_ABS))
        else
            push!(diss_scores, clip01(1 - float(Dex) / D_REF_ABS))
        end
        δ = get(c, "shock_thickness", nothing)
        if δ === nothing || δ === missing || !isfinite(float(δ))
            sδ = 1.0
        else
            sδ = clip01(1 - float(δ) / DELTA_REF_ABS)
        end
        sη = clip01(1 - η / ETA_REF_ABS)
        push!(shock_scores, 0.5 * sδ + 0.5 * sη)
    end
    S_diss = isempty(diss_scores) ? 1.0 : sum(diss_scores) / length(diss_scores)
    S_shock = isempty(shock_scores) ? 1.0 : sum(shock_scores) / length(shock_scores)

    # Robustness
    robust_fail = any(
        c ->
            get(c, "diverged", false) === true ||
            get(c, "nan_detected", false) === true ||
            get(c, "positivity_ok", true) === false,
        cases,
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
    apply_scores!(report; baseline=nothing)

Fill `report["summary"]["scores"]` using absolute maps.
(Relative maps are applied in invent classification, not here.)
"""
function apply_scores!(report::AbstractDict; baseline=nothing)
    cases = get(report, "cases", Any[])
    scores = score_suite_absolute(cases)
    report["summary"]["scores"] = scores
    report["scoring_formula_version"] = SCORING_FORMULA_VERSION
    return report
end

"""Hard-gate failures from cases (for suite orchestration)."""
function collect_hard_gate_failures(cases::AbstractVector)
    fails = String[]
    for c in cases
        name = get(c, "name", "?")
        if get(c, "diverged", false) || get(c, "nan_detected", false)
            push!(fails, "$name diverged/nan")
        end
        if get(c, "case_type", "") == "smooth_order" && !get(c, "order_pass", true)
            push!(fails, "$name order_pass failed")
        end
        if get(c, "equation", "") == "euler1d" && get(c, "positivity_ok", true) === false
            push!(fails, "$name positivity failed")
        end
        if get(c, "conservation_metric", "none") != "none" &&
           !get(c, "conservation_pass", true)
            push!(fails, "$name conservation_pass failed")
        end
    end
    return fails
end
