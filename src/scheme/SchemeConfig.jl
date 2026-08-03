# Base-scheme configuration (solution points, numerical flux, time integrator).
#
# Invent loop freezes DEFAULT_SCHEME / FROZEN_INVENT_SCHEME (GL + Rusanov + SSP-RK3)
# so composite scores remain comparable across history unless a logged re-baseline.

"""
    SchemeConfig

Discrete scheme axes orthogonal to the capturing method.

- `points`: `:gl` (default) or `:gll`
- `flux`: `:rusanov` (default) or `:hllc`
- `time`: `:ssp_rk3` (default) or `:ssp_rk2`
"""
struct SchemeConfig
    points::Symbol
    flux::Symbol
    time::Symbol
end

function SchemeConfig(;
    points::Symbol=:gl,
    flux::Symbol=:rusanov,
    time::Symbol=:ssp_rk3,
)
    points in (:gl, :gll) || throw(ArgumentError("points must be :gl or :gll, got $points"))
    flux in (:rusanov, :hllc) || throw(ArgumentError("flux must be :rusanov or :hllc, got $flux"))
    time in (:ssp_rk3, :ssp_rk2) || throw(ArgumentError("time must be :ssp_rk3 or :ssp_rk2, got $time"))
    return SchemeConfig(points, flux, time)
end

"""Default / invent-frozen scheme: GL + Rusanov + SSP-RK3."""
const DEFAULT_SCHEME = SchemeConfig(:gl, :rusanov, :ssp_rk3)

"""
    scheme_dict(scheme) -> Dict{String,Any}

JSON-friendly scheme record (human labels for reports / experiment log).
"""
function scheme_dict(scheme::SchemeConfig)
    return Dict{String,Any}(
        "points" => scheme.points === :gl ? "GL" : "GLL",
        "flux" => scheme.flux === :rusanov ? "Rusanov" : "HLLC",
        "time" => scheme.time === :ssp_rk3 ? "SSP-RK3" : "SSP-RK2",
        "points_symbol" => string(scheme.points),
        "flux_symbol" => string(scheme.flux),
        "time_symbol" => string(scheme.time),
    )
end

"""
    parse_scheme(; points="gl", flux="rusanov", time="ssp_rk3") -> SchemeConfig

Parse CLI / string scheme axes (case-insensitive).
"""
function parse_scheme(;
    points::AbstractString="gl",
    flux::AbstractString="rusanov",
    time::AbstractString="ssp_rk3",
)
    p = Symbol(lowercase(String(points)))
    f = Symbol(lowercase(String(flux)))
    t = Symbol(lowercase(replace(String(time), "-" => "_")))
    return SchemeConfig(; points=p, flux=f, time=t)
end

"""CFL guidance note for time integrators (documentation / reports)."""
function time_cfl_guidance(time::Symbol)
    if time === :ssp_rk3
        return "SSP-RK3: default CFL ≈ 0.2 (order studies may shrink CFL with p)"
    elseif time === :ssp_rk2
        return "SSP-RK2: same CFL form as SSP-RK3; typically use ≤ SSP-RK3 CFL (default 0.2); second-order in time — prefer fixed small Δt for high-order spatial studies"
    else
        return "unknown time scheme"
    end
end
