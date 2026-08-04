# 2D compressible Euler (4 conserved variables).

"""
    Euler2D{T}(γ=1.4)

Conserved U = (ρ, ρu, ρv, E). Reuses 1D EOS helpers where possible.
"""
struct Euler2D{T} <: AbstractEquation{4}
    γ::T
end

Euler2D(γ::Real=1.4) = Euler2D{typeof(float(γ))}(float(γ))

function pressure(eq::Euler2D{T}, U::AbstractVector) where {T}
    ρ = max(U[1], eps(T))
    u = U[2] / ρ
    v = U[3] / ρ
    return (eq.γ - one(T)) * (U[4] - T(0.5) * ρ * (u * u + v * v))
end

function primitives_to_conserved(eq::Euler2D{T}, ρ, u, v, p) where {T}
    ρT, uT, vT, pT = T(ρ), T(u), T(v), T(p)
    E = pT / (eq.γ - one(T)) + T(0.5) * ρT * (uT * uT + vT * vT)
    return T[ρT, ρT * uT, ρT * vT, E]
end

function conserved_to_primitives(eq::Euler2D{T}, U::AbstractVector) where {T}
    ρ = U[1]
    u = U[2] / max(ρ, eps(T))
    v = U[3] / max(ρ, eps(T))
    p = pressure(eq, U)
    return ρ, u, v, p
end

function sound_speed(eq::Euler2D{T}, U::AbstractVector) where {T}
    ρ = max(U[1], eps(T))
    p = max(pressure(eq, U), eps(T))
    return sqrt(eq.γ * p / ρ)
end

function physical_flux_x(eq::Euler2D{T}, U::AbstractVector) where {T}
    out = Vector{T}(undef, 4)
    physical_flux_x!(out, eq, U)
    return out
end

function physical_flux_y(eq::Euler2D{T}, U::AbstractVector) where {T}
    out = Vector{T}(undef, 4)
    physical_flux_y!(out, eq, U)
    return out
end

"""In-place physical flux (x). `out` length ≥ 4."""
function physical_flux_x!(out::AbstractVector{T}, eq::Euler2D{T}, U::AbstractVector) where {T}
    ρ, ρu, ρv, E = U[1], U[2], U[3], U[4]
    u = ρu / max(ρ, eps(T))
    p = pressure(eq, U)
    @inbounds begin
        out[1] = ρu
        out[2] = ρu * u + p
        out[3] = ρv * u
        out[4] = (E + p) * u
    end
    return out
end

"""In-place physical flux (y). `out` length ≥ 4."""
function physical_flux_y!(out::AbstractVector{T}, eq::Euler2D{T}, U::AbstractVector) where {T}
    ρ, ρu, ρv, E = U[1], U[2], U[3], U[4]
    v = ρv / max(ρ, eps(T))
    p = pressure(eq, U)
    @inbounds begin
        out[1] = ρv
        out[2] = ρu * v
        out[3] = ρv * v + p
        out[4] = (E + p) * v
    end
    return out
end

function max_wave_speed_n(eq::Euler2D, uL::AbstractVector, uR::AbstractVector, nx, ny)
    function spd(U)
        ρ = max(U[1], eps(eltype(U)))
        u = U[2] / ρ
        v = U[3] / ρ
        un = u * nx + v * ny
        return abs(un) + sound_speed(eq, U)
    end
    return max(spd(uL), spd(uR))
end

function max_wave_speed(eq::Euler2D{T}, U::AbstractArray{T,4}) where {T}
    m = zero(T)
    Np1, Np2, Nel = size(U, 1), size(U, 2), size(U, 3)
    @inbounds for e in 1:Nel, j in 1:Np2, i in 1:Np1
        Uv = @view U[i, j, e, :]
        ρ = max(Uv[1], eps(T))
        u = Uv[2] / ρ
        v = Uv[3] / ρ
        c = sound_speed(eq, Uv)
        m = max(m, abs(u) + c, abs(v) + c)
    end
    return max(m, eps(T))
end

"""Rusanov numerical flux for normal direction n=(nx,ny)."""
function numerical_flux_n(eq::Euler2D{T}, uL::AbstractVector, uR::AbstractVector, nx, ny) where {T}
    out = Vector{T}(undef, 4)
    numerical_flux_n!(out, eq, uL, uR, nx, ny)
    return out
end

"""In-place Rusanov normal flux. Allocates nothing beyond temporaries of size Neq."""
function numerical_flux_n!(
    out::AbstractVector{T},
    eq::Euler2D{T},
    uL::AbstractVector,
    uR::AbstractVector,
    nx,
    ny,
) where {T}
    # F·n = nx Fx + ny Gy
    ρL = max(uL[1], eps(T))
    uL_ = uL[2] / ρL
    vL_ = uL[3] / ρL
    pL = pressure(eq, uL)
    ρR = max(uR[1], eps(T))
    uR_ = uR[2] / ρR
    vR_ = uR[3] / ρR
    pR = pressure(eq, uR)
    nxT, nyT = T(nx), T(ny)
    # Physical flux · n
    fL1 = uL[2] * nxT + uL[3] * nyT
    fL2 = (uL[2] * uL_ + pL) * nxT + (uL[2] * vL_) * nyT
    fL3 = (uL[3] * uL_) * nxT + (uL[3] * vL_ + pL) * nyT
    fL4 = (uL[4] + pL) * (uL_ * nxT + vL_ * nyT)
    fR1 = uR[2] * nxT + uR[3] * nyT
    fR2 = (uR[2] * uR_ + pR) * nxT + (uR[2] * vR_) * nyT
    fR3 = (uR[3] * uR_) * nxT + (uR[3] * vR_ + pR) * nyT
    fR4 = (uR[4] + pR) * (uR_ * nxT + vR_ * nyT)
    λ = max_wave_speed_n(eq, uL, uR, nx, ny)
    half = T(0.5)
    @inbounds begin
        out[1] = half * (fL1 + fR1) - half * T(λ) * (uR[1] - uL[1])
        out[2] = half * (fL2 + fR2) - half * T(λ) * (uR[2] - uL[2])
        out[3] = half * (fL3 + fR3) - half * T(λ) * (uR[3] - uL[3])
        out[4] = half * (fL4 + fR4) - half * T(λ) * (uR[4] - uL[4])
    end
    return out
end

# Allow positivity_ok(eq::Euler2D, state) via Euler1D-style on 4-component
function positivity_ok(eq::Euler2D{T}, state::SolutionState2D{T,4}; atol=zero(T)) where {T}
    Np = size(state.u, 1)
    for e in 1:state.mesh.n_elements, j in 1:Np, i in 1:Np
        U = @view state.u[i, j, e, :]
        p = pressure(eq, U)
        if !(U[1] > atol && p > atol && isfinite(U[1]) && isfinite(p))
            return false
        end
    end
    return true
end
