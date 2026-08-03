# Order and conservation metrics for verification.

const DEFAULT_ORDER_TOLERANCE = 0.3
const DEFAULT_CONSERVATION_TOL_REL = 1e-12

"""
    observed_orders(mesh_sizes, errors) -> Vector

For successive refinements, q_i = log2(E_i / E_{i+1}) when h halves.
`mesh_sizes` are characteristic h (e.g. Δx or 1/N_el); must be refining.
"""
function observed_orders(mesh_sizes::AbstractVector, errors::AbstractVector)
    n = length(errors)
    n >= 2 || return Float64[]
    qs = Float64[]
    for i in 1:(n - 1)
        h0, h1 = mesh_sizes[i], mesh_sizes[i + 1]
        e0, e1 = errors[i], errors[i + 1]
        # Prefer log2(E0/E1) when h is halved; general: log(E0/E1)/log(h0/h1)
        if e0 <= 0 || e1 <= 0
            push!(qs, NaN)
        else
            push!(qs, log(e0 / e1) / log(h0 / h1))
        end
    end
    return qs
end

"""
    order_pass(observed, formal_order; tol=0.3) -> Bool

True if every finite observed order q satisfies q >= formal_order - tol.
"""
function order_pass(
    observed::AbstractVector;
    formal_order::Real,
    tol::Real=DEFAULT_ORDER_TOLERANCE,
)
    isempty(observed) && return false
    for q in observed
        !isfinite(q) && return false
        q < formal_order - tol && return false
    end
    return true
end

"""
Relative conservation residual |M(T)-M(0)| / scale.

When the conserved integral is near zero (e.g. pure sine on a period),
relative residual is ill-defined; we scale by `max(|M0|, |MT|, 1)` so the
metric reports an absolute change of O(ε) as ~ε.
"""
function conservation_residual_relative(M0, MT)
    T = typeof(float(M0))
    denom = max(abs(M0), abs(MT), one(T))
    return abs(MT - M0) / denom
end

"""Absolute mass change |M(T)-M(0)|."""
conservation_residual_absolute(M0, MT) = abs(MT - M0)

"""Min/max of solution field (component `c`)."""
function solution_extrema(state::SolutionState, c::Int=1)
    umin = typemax(eltype(state.u))
    umax = typemin(eltype(state.u))
    @inbounds for e in 1:size(state.u, 2), j in 1:size(state.u, 1)
        v = state.u[j, e, c]
        umin = min(umin, v)
        umax = max(umax, v)
    end
    return umin, umax
end

"""
    overshoot_metric(u_min, u_max, u_ref_min, u_ref_max) -> (η_over, η_under, η)

Normalized overshoot / undershoot relative to reference bounds (e.g. IC range).
Jump scale = max(|u_ref_max - u_ref_min|, eps).
"""
function overshoot_metric(u_min, u_max, u_ref_min, u_ref_max)
    Δ = max(abs(u_ref_max - u_ref_min), eps(Float64))
    η_over = max(0.0, float(u_max) - float(u_ref_max)) / Δ
    η_under = max(0.0, float(u_ref_min) - float(u_min)) / Δ
    η = max(η_over, η_under)
    return η_over, η_under, η
end
