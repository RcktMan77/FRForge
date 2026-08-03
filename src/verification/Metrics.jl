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

"""
    sample_solution_1d(state, eq; component=:density) -> (x, v)

Collect all solution-point samples sorted by physical x.
`component` is `:density`, `:velocity`, `:pressure`, or an Int conserved index.
"""
function sample_solution_1d(state::SolutionState{T}, eq; component=:density) where {T}
    mesh, ops = state.mesh, state.ops
    Np, Nel = size(state.u, 1), size(state.u, 2)
    n = Np * Nel
    xs = Vector{T}(undef, n)
    vs = Vector{T}(undef, n)
    idx = 1
    for e in 1:Nel
        xe = physical_coords(mesh, ops, e)
        for j in 1:Np
            U = @view state.u[j, e, :]
            xs[idx] = xe[j]
            if component === :density || component === 1
                vs[idx] = U[1]
            elseif component === :velocity
                vs[idx] = velocity(eq, U)
            elseif component === :pressure
                vs[idx] = pressure(eq, U)
            elseif component isa Integer
                vs[idx] = U[component]
            else
                error("Unknown component $component")
            end
            idx += 1
        end
    end
    perm = sortperm(xs)
    return xs[perm], vs[perm]
end

"""Mean solution-point spacing along the sorted chain."""
function mean_sp_spacing(x::AbstractVector)
    n = length(x)
    n < 2 && return 1.0
    return (x[end] - x[1]) / (n - 1)
end

"""
    shock_thickness_sp(x, ρ; jump_lo=nothing, jump_hi=nothing) -> δ

10%–90% rise distance of the primary density jump, in units of mean SP spacing.
If jump bounds are not given, uses global min/max of ρ.
"""
function shock_thickness_sp(
    x::AbstractVector,
    ρ::AbstractVector;
    jump_lo=nothing,
    jump_hi=nothing,
    window_sp::Int=30,
)
    n = length(x)
    n < 4 && return NaN
    # Find steepest gradient location (primary shock locus)
    imax = 2
    gmax = 0.0
    for i in 2:n
        g = abs(ρ[i] - ρ[i - 1]) / max(abs(x[i] - x[i - 1]), 1e-30)
        if g > gmax
            gmax = g
            imax = i
        end
    end
    # Local window around shock to estimate left/right plateaus
    i0 = max(1, imax - window_sp)
    i1 = min(n, imax + window_sp)
    if jump_lo === nothing || jump_hi === nothing
        # median of left third and right third of window
        w = i1 - i0 + 1
        nL = max(1, w ÷ 3)
        left_vals = sort(ρ[i0:(i0 + nL - 1)])
        right_vals = sort(ρ[(i1 - nL + 1):i1])
        ρL = left_vals[max(1, length(left_vals) ÷ 2)]
        ρR = right_vals[max(1, length(right_vals) ÷ 2)]
        ρmin, ρmax = extrema((ρL, ρR))
    else
        ρmin = float(jump_lo)
        ρmax = float(jump_hi)
    end
    Δ = ρmax - ρmin
    abs(Δ) < 1e-14 && return NaN
    lo = ρmin + 0.1 * Δ
    hi = ρmin + 0.9 * Δ
    rising = ρR >= ρL
    i10, i90 = imax, imax
    if rising
        i10 = imax
        while i10 > i0 && ρ[i10] > lo
            i10 -= 1
        end
        i90 = imax
        while i90 < i1 && ρ[i90] < hi
            i90 += 1
        end
    else
        i10 = imax
        while i10 < i1 && ρ[i10] > lo
            i10 += 1
        end
        i90 = imax
        while i90 > i0 && ρ[i90] < hi
            i90 -= 1
        end
    end
    dx = abs(x[i90] - x[i10])
    return dx / mean_sp_spacing(x)
end
"""
    excess_dissipation(x, u_method, u_ref; mask=nothing) -> D_ex

Relative L1 difference on mask S:
  D_ex = Σ_{S} |u_m - u_ref| / (Σ_S |u_ref| + ε)
If `mask` is nothing, use all points.
"""
function excess_dissipation(
    x::AbstractVector,
    u_method::AbstractVector,
    u_ref::AbstractVector;
    mask::Union{Nothing,AbstractVector{Bool}}=nothing,
)
    n = length(u_method)
    num = 0.0
    den = 0.0
    for i in 1:n
        (mask !== nothing && !mask[i]) && continue
        num += abs(u_method[i] - u_ref[i])
        den += abs(u_ref[i])
    end
    return num / (den + 1e-30)
end

"""
Smooth-region mask: |∂ρ/∂x| below fraction of max gradient, optionally excluding windows.
"""
function smooth_region_mask(
    x::AbstractVector,
    ρ::AbstractVector;
    grad_frac::Float64=0.1,
    exclude_windows::Vector{Tuple{Float64,Float64}}=Tuple{Float64,Float64}[],
)
    n = length(x)
    mask = trues(n)
    # gradients
    g = zeros(n)
    for i in 2:n
        g[i] = abs(ρ[i] - ρ[i - 1]) / max(abs(x[i] - x[i - 1]), 1e-30)
    end
    g[1] = g[2]
    gmax = maximum(g)
    thr = grad_frac * gmax
    for i in 1:n
        if g[i] > thr
            mask[i] = false
        end
        for (a, b) in exclude_windows
            if a <= x[i] <= b
                mask[i] = false
            end
        end
    end
    return mask
end

"""Linear interpolate reference (x_ref, v_ref) onto query points xq."""
function interp_1d(x_ref::AbstractVector, v_ref::AbstractVector, xq::AbstractVector)
    n = length(xq)
    out = similar(xq, Float64)
    for i in 1:n
        x = float(xq[i])
        if x <= x_ref[1]
            out[i] = float(v_ref[1])
        elseif x >= x_ref[end]
            out[i] = float(v_ref[end])
        else
            j = searchsortedlast(x_ref, x)
            j = clamp(j, 1, length(x_ref) - 1)
            t = (x - x_ref[j]) / (x_ref[j + 1] - x_ref[j] + 1e-30)
            out[i] = (1 - t) * float(v_ref[j]) + t * float(v_ref[j + 1])
        end
    end
    return out
end
