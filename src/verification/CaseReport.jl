# Shared verification case-dict fields (schema v1 keys). Behavior-preserving helpers only.

"""
    case_report_dict(; name, case_type, equation, p, capturing_method, pass,
                       diverged, nan_detected, kwargs...) -> Dict{String,Any}

Build a case report with the standard required keys, plus optional extra fields.
"""
function case_report_dict(;
    name::AbstractString,
    case_type::AbstractString,
    equation::AbstractString,
    p::Integer,
    capturing_method::AbstractString="null",
    pass::Bool,
    diverged::Bool=false,
    nan_detected::Bool=false,
    conservation_residual=0.0,
    conservation_pass::Bool=true,
    conservation_metric::AbstractString="none",
    positivity_ok::Bool=true,
    wall_time_sec::Real=0.0,
    n_elements::Union{Nothing,Integer}=nothing,
    t_final=nothing,
    excess_dissipation=nothing,
    shock_thickness=nothing,
    shock_thickness_unit::AbstractString="sp_spacings",
    overshoot=0.0,
    metrics::Union{Nothing,AbstractDict}=nothing,
    extra::AbstractDict=Dict{String,Any}(),
)
    d = Dict{String,Any}(
        "name" => String(name),
        "case_type" => String(case_type),
        "equation" => String(equation),
        "p" => Int(p),
        "capturing_method" => String(capturing_method),
        "pass" => pass,
        "diverged" => diverged,
        "nan_detected" => nan_detected,
        "conservation_residual" => conservation_residual,
        "conservation_pass" => conservation_pass,
        "conservation_metric" => String(conservation_metric),
        "positivity_ok" => positivity_ok,
        "wall_time_sec" => Float64(wall_time_sec),
        "excess_dissipation" => excess_dissipation,
        "shock_thickness" => shock_thickness,
        "shock_thickness_unit" => String(shock_thickness_unit),
        "overshoot" => overshoot,
    )
    n_elements !== nothing && (d["n_elements"] = Int(n_elements))
    t_final !== nothing && (d["t_final"] = t_final)
    metrics !== nothing && (d["metrics"] = Dict{String,Any}(metrics))
    for (k, v) in pairs(extra)
        d[String(k)] = v
    end
    return d
end
