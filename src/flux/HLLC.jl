# HLLC approximate Riemann solver for compressible Euler (Toro).
# Used as the less-dissipative flux corner in the Phase 2 robustness matrix.

"""
    hllc_flux(eq::Euler1D, uL, uR) -> numerical flux

1D HLLC with direct (Davis) wave-speed estimates.
"""
function hllc_flux(eq::Euler1D{T}, uL::AbstractVector, uR::AbstractVector) where {T}
    γ = eq.γ
    ρL = max(T(uL[1]), eps(T))
    ρR = max(T(uR[1]), eps(T))
    uvel_L = T(uL[2]) / ρL
    uvel_R = T(uR[2]) / ρR
    pL = max(pressure(eq, uL), eps(T))
    pR = max(pressure(eq, uR), eps(T))
    cL = sqrt(γ * pL / ρL)
    cR = sqrt(γ * pR / ρR)

    SL = min(uvel_L - cL, uvel_R - cR)
    SR = max(uvel_L + cL, uvel_R + cR)

    # Contact speed S_M
    num = pR - pL + ρL * uvel_L * (SL - uvel_L) - ρR * uvel_R * (SR - uvel_R)
    den = ρL * (SL - uvel_L) - ρR * (SR - uvel_R)
    if abs(den) < eps(T)
        SM = T(0.5) * (uvel_L + uvel_R)
    else
        SM = num / den
    end

    FL = physical_flux(eq, uL)
    FR = physical_flux(eq, uR)

    if zero(T) <= SL
        return collect(T, FL)
    elseif zero(T) >= SR
        return collect(T, FR)
    end

    function star_state(U, ρ, uvel, p, S)
        # U* = ρ(S-u)/(S-S_M) * [1, S_M, E/ρ + (S_M-u)(S_M + p/(ρ(S-u)))]
        EL = U[3]
        factor = ρ * (S - uvel) / (S - SM)
        e_int = EL / ρ + (SM - uvel) * (SM + p / (ρ * (S - uvel)))
        return T[factor, factor * SM, factor * e_int]
    end

    if SL <= zero(T) <= SM
        Ustar = star_state(uL, ρL, uvel_L, pL, SL)
        out = Vector{T}(undef, 3)
        @inbounds for c in 1:3
            out[c] = FL[c] + SL * (Ustar[c] - uL[c])
        end
        return out
    else
        # SM <= 0 <= SR
        Ustar = star_state(uR, ρR, uvel_R, pR, SR)
        out = Vector{T}(undef, 3)
        @inbounds for c in 1:3
            out[c] = FR[c] + SR * (Ustar[c] - uR[c])
        end
        return out
    end
end

"""
    hllc_flux_n(eq::Euler2D, uL, uR, nx, ny) -> normal numerical flux

HLLC in the face-normal direction for 2D Euler.
"""
function hllc_flux_n(
    eq::Euler2D{T},
    uL::AbstractVector,
    uR::AbstractVector,
    nx,
    ny,
) where {T}
    nxT, nyT = T(nx), T(ny)
    # Normalize (safety)
    nrm = hypot(nxT, nyT)
    nrm < eps(T) && return numerical_flux_n(eq, uL, uR, nx, ny)
    nxT /= nrm
    nyT /= nrm

    γ = eq.γ
    function prim(U)
        ρ = max(U[1], eps(T))
        u = U[2] / ρ
        v = U[3] / ρ
        p = max(pressure(eq, U), eps(T))
        un = u * nxT + v * nyT
        return ρ, u, v, p, un
    end
    ρL, uL_, vL_, pL, unL = prim(uL)
    ρR, uR_, vR_, pR, unR = prim(uR)
    cL = sqrt(γ * pL / ρL)
    cR = sqrt(γ * pR / ρR)

    SL = min(unL - cL, unR - cR)
    SR = max(unL + cL, unR + cR)
    num = pR - pL + ρL * unL * (SL - unL) - ρR * unR * (SR - unR)
    den = ρL * (SL - unL) - ρR * (SR - unR)
    SM = abs(den) < eps(T) ? T(0.5) * (unL + unR) : num / den

    fL = physical_flux_x(eq, uL) .* nxT .+ physical_flux_y(eq, uL) .* nyT
    fR = physical_flux_x(eq, uR) .* nxT .+ physical_flux_y(eq, uR) .* nyT

    if zero(T) <= SL
        return collect(T, fL)
    elseif zero(T) >= SR
        return collect(T, fR)
    end

    function star_state(U, ρ, u, v, p, un, S)
        # Tangent velocity unchanged; normal velocity → SM
        ut_x = u - un * nxT
        ut_y = v - un * nyT
        u_star = ut_x + SM * nxT
        v_star = ut_y + SM * nyT
        factor = ρ * (S - un) / (S - SM)
        E = U[4]
        e_star = E / ρ + (SM - un) * (SM + p / (ρ * (S - un)))
        return T[factor, factor * u_star, factor * v_star, factor * e_star]
    end

    if SL <= zero(T) <= SM
        Us = star_state(uL, ρL, uL_, vL_, pL, unL, SL)
        out = Vector{T}(undef, 4)
        @inbounds for c in 1:4
            out[c] = fL[c] + SL * (Us[c] - uL[c])
        end
        return out
    else
        Us = star_state(uR, ρR, uR_, vR_, pR, unR, SR)
        out = Vector{T}(undef, 4)
        @inbounds for c in 1:4
            out[c] = fR[c] + SR * (Us[c] - uR[c])
        end
        return out
    end
end

"""
    interface_flux(eq, uL, uR, flux_kind) -> flux

Dispatch numerical flux by scheme axis. `:rusanov` uses each equation's
default `numerical_flux` (preserves pure upwind for linear advection).
`:hllc` uses HLLC for Euler; falls back to equation default otherwise.
"""
function interface_flux(eq::AbstractEquation, uL, uR, flux_kind::Symbol)
    if flux_kind === :hllc
        return _hllc_or_fallback(eq, uL, uR)
    elseif flux_kind === :rusanov
        return numerical_flux(eq, uL, uR)
    else
        throw(ArgumentError("Unknown flux kind $flux_kind"))
    end
end

function _hllc_or_fallback(eq::Euler1D, uL::AbstractVector, uR::AbstractVector)
    return hllc_flux(eq, uL, uR)
end

function _hllc_or_fallback(eq::AbstractEquation, uL, uR)
    # Scalar laws / non-Euler: HLLC not defined — use equation default
    return numerical_flux(eq, uL, uR)
end

"""Normal-interface flux for 2D (scheme-aware)."""
function interface_flux_n(eq::AbstractEquation, uL, uR, nx, ny, flux_kind::Symbol)
    if flux_kind === :hllc
        return _hllc_or_fallback_n(eq, uL, uR, nx, ny)
    elseif flux_kind === :rusanov
        return numerical_flux_n(eq, uL, uR, nx, ny)
    else
        throw(ArgumentError("Unknown flux kind $flux_kind"))
    end
end

"""In-place normal-interface flux (Rusanov path; HLLC falls back to allocating)."""
function interface_flux_n!(out::AbstractVector, eq::AbstractEquation, uL, uR, nx, ny, flux_kind::Symbol)
    if flux_kind === :rusanov
        return numerical_flux_n!(out, eq, uL, uR, nx, ny)
    else
        fh = interface_flux_n(eq, uL, uR, nx, ny, flux_kind)
        @inbounds for c in 1:length(out)
            out[c] = fh[c]
        end
        return out
    end
end

function _hllc_or_fallback_n(eq::Euler2D, uL, uR, nx, ny)
    return hllc_flux_n(eq, uL, uR, nx, ny)
end

function _hllc_or_fallback_n(eq::AbstractEquation, uL, uR, nx, ny)
    return numerical_flux_n(eq, uL, uR, nx, ny)
end
