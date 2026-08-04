# Shared pure helpers for Persson-style sensors and BR0 AV (1D/2D).
# Extracted for maintainability; formulas match existing PerssonAV / PerssonAV2D.

"""
    persson_s0(p, s0_factor) -> s0

Threshold for the smooth modal indicator: `s0_factor * log10(max(p,1))`.
"""
@inline function persson_s0(p::Integer, s0_factor::T) where {T<:Real}
    return s0_factor * log10(T(max(p, 1)))
end

"""
    persson_sigma(e_high, e_tot, s0, κ, ε_floor) -> σ ∈ (0,1)

Smooth logistic indicator from modal energies (same expression as before).
"""
@inline function persson_sigma(e_high::T, e_tot::T, s0::T, κ::T, ε_floor::T) where {T<:Real}
    s_e = log10(e_high / e_tot + ε_floor)
    return one(T) / (one(T) + exp(-κ * (s_e - s0)))
end

"""
    br0_tau_1d(εm, εp, hm, hp, p) -> τ

1D BR0 interior penalty (moderate form used by `apply_dissipation_br0!`):
`τ = 2 * ε̄ * (p+1) / h̄`.
"""
@inline function br0_tau_1d(εm::T, εp::T, hm::T, hp::T, p::Integer) where {T<:Real}
    εbar = T(0.5) * (εm + εp)
    hbar = T(0.5) * (hm + hp)
    return T(2) * εbar * T(p + 1) / max(hbar, eps(T))
end

"""
    br0_avg_flux(gm, gp, um, up, τ) -> ĝ

Central average minus penalty: `½(gm+gp) − τ (up−um)`.
"""
@inline function br0_avg_flux(gm::T, gp::T, um::T, up::T, τ::T) where {T<:Real}
    return T(0.5) * (gm + gp) - τ * (up - um)
end
