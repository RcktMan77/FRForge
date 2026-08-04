# Coefficient-tuning helper (secondary invent path).
#
# Structural invent against residual hooks remains primary. This module only
# grids a scalar coefficient of a *registered factory* on the light quant suite
# for quick scouting — it does **not** change DEFAULT_SCHEME or invent history
# unless the caller explicitly invents a chosen point.

"""
    tune_coefficient(base_method; param=:c_av, values, baseline="persson_av",
                     suite=:light, δ=DEFAULT_SCORE_MARGIN) -> Vector{Dict}

Run the light (or `:quant`) suite for each coefficient in `values` by registering
a temporary method that wraps `base_method` factory kwargs.

Returns one row per value with composite scores and candidate_status vs baseline.
Does not append the experiment log (caller may invent the best point later).
"""
function tune_coefficient(
    base_method::AbstractString;
    param::Symbol=:c_av,
    values::AbstractVector{<:Real}=[0.05, 0.1, 0.2, 0.5],
    baseline::AbstractString="persson_av",
    suite::Symbol=:light,
    δ::Real=DEFAULT_SCORE_MARGIN,
    seed=nothing,
)
    base_method = require_registered_method(base_method)
    baseline = require_registered_method(baseline; role="baseline")
    haskey(METHOD_REGISTRY, base_method) || error("unknown base method $base_method")

    bas = run_method_report(baseline; suite=suite, seed=seed)
    rows = Dict{String,Any}[]
    println("tune: base=$base_method param=$param values=$values suite=$suite baseline=$baseline")
    println("  (does not append experiment log; invent frozen scheme unchanged)")

    for (i, v) in enumerate(values)
        tmp_name = string("_tune_", base_method, "_", param, "_", i)
        # Factory: call registered base factory with param override when supported
        register_method!(
            tmp_name,
            (; kwargs...) -> begin
                kw = Dict{Symbol,Any}(pairs(kwargs))
                kw[param] = v
                return METHOD_REGISTRY[base_method](; kw...)
            end,
        )
        try
            met = run_method_report(tmp_name; suite=suite, seed=seed)
            cmp = classify_candidate(met, bas; δ=δ, vtk_produced=false)
            row = Dict{String,Any}(
                "param" => string(param),
                "value" => float(v),
                "tmp_method" => tmp_name,
                "candidate_status" => cmp["candidate_status"],
                "composite" => cmp["absolute_scores"]["composite"],
                "baseline_composite" => cmp["baseline_composite"],
                "composite_margin" => cmp["composite_margin"],
                "tradeoff_ok" => cmp["tradeoff_ok"],
                "overall_pass" => get(met, "overall_pass", false),
            )
            push!(rows, row)
            println(
                "  $param=$v  status=$(row["candidate_status"])  composite=$(row["composite"])  margin=$(row["composite_margin"])",
            )
        finally
            delete!(METHOD_REGISTRY, tmp_name)
        end
    end

    # rank by composite descending among overall_pass
    sort!(rows; by=r -> (r["overall_pass"] === true ? 0 : 1, -float(r["composite"])))
    if !isempty(rows)
        best = rows[1]
        println(
            "best scout: $(best["param"])=$(best["value"]) status=$(best["candidate_status"]) composite=$(best["composite"])",
        )
        println("  next: register a method with that coefficient, then `frforge invent` for log history")
    end
    return rows
end
