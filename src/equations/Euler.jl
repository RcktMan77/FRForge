# 1D compressible Euler equations (conservative form).

"""
    Euler1D{T}(γ=1.4)

1D Euler: conserved variables U = (ρ, ρu, E).
Ideal gas with ratio of specific heats `γ`.
"""
struct Euler1D{T} <: AbstractEquation{3}
    γ::T
end

Euler1D(γ::Real=1.4) = Euler1D{typeof(float(γ))}(float(γ))
Euler1D{T}() where {T} = Euler1D{T}(T(1.4))

# --- Primitive / conserved conversions ---

"""Pressure from conserved state (ρ, ρu, E)."""
function pressure(eq::Euler1D{T}, U::AbstractVector) where {T}
    ρ, ρu, E = U[1], U[2], U[3]
    u = ρu / ρ
    return (eq.γ - one(T)) * (E - T(0.5) * ρ * u * u)
end

function pressure(eq::Euler1D{T}, ρ, ρu, E) where {T}
    u = ρu / ρ
    return (eq.γ - one(T)) * (E - T(0.5) * ρ * u * u)
end

"""Velocity from conserved state."""
function velocity(::Euler1D, U::AbstractVector)
    return U[2] / U[1]
end

"""Sound speed c = √(γ p / ρ)."""
function sound_speed(eq::Euler1D{T}, U::AbstractVector) where {T}
    ρ = U[1]
    p = pressure(eq, U)
    return sqrt(eq.γ * p / ρ)
end

"""Build conserved state from primitives (ρ, u, p)."""
function primitives_to_conserved(eq::Euler1D{T}, ρ, u, p) where {T}
    ρT = T(ρ)
    uT = T(u)
    pT = T(p)
    E = pT / (eq.γ - one(T)) + T(0.5) * ρT * uT * uT
    return T[ρT, ρT * uT, E]
end

"""Primitives (ρ, u, p) from conserved U."""
function conserved_to_primitives(eq::Euler1D{T}, U::AbstractVector) where {T}
    ρ = U[1]
    u = U[2] / ρ
    p = pressure(eq, U)
    return ρ, u, p
end

"""Positivity: ρ > 0 and p > 0 for a single conserved state."""
function positivity_ok_state(eq::Euler1D{T}, U::AbstractVector; atol=zero(T)) where {T}
    ρ = U[1]
    p = pressure(eq, U)
    return ρ > atol && p > atol && isfinite(ρ) && isfinite(p)
end

"""Positivity over a 3D solution array u[j,e,c]."""
function positivity_ok(eq::Euler1D{T}, u::AbstractArray{T,3}; atol=zero(T)) where {T}
    Np, Nel = size(u, 1), size(u, 2)
    @inbounds for e in 1:Nel, j in 1:Np
        if !positivity_ok_state(eq, @view(u[j, e, :]); atol=atol)
            return false
        end
    end
    return true
end

# --- Flux and wave speeds ---

function physical_flux(eq::Euler1D{T}, U::AbstractVector) where {T}
    ρ, ρu, E = U[1], U[2], U[3]
    u = ρu / ρ
    p = pressure(eq, U)
    return T[ρu, ρu * u + p, (E + p) * u]
end

function physical_flux(eq::Euler1D{T}, Umat::AbstractMatrix) where {T}
    Np = size(Umat, 1)
    F = similar(Umat)
    @inbounds for j in 1:Np
        f = physical_flux(eq, @view Umat[j, :])
        F[j, 1] = f[1]
        F[j, 2] = f[2]
        F[j, 3] = f[3]
    end
    return F
end

function max_wave_speed(eq::Euler1D, uL::AbstractVector, uR::AbstractVector)
    vL = velocity(eq, uL)
    vR = velocity(eq, uR)
    cL = sound_speed(eq, uL)
    cR = sound_speed(eq, uR)
    return max(abs(vL) + cL, abs(vR) + cR)
end

function max_wave_speed(eq::Euler1D{T}, U::AbstractArray{T,3}) where {T}
    m = zero(T)
    Np, Nel = size(U, 1), size(U, 2)
    @inbounds for e in 1:Nel, j in 1:Np
        Uj = @view U[j, e, :]
        v = abs(velocity(eq, Uj)) + sound_speed(eq, Uj)
        m = max(m, v)
    end
    return max(m, eps(T))
end

function numerical_flux(eq::Euler1D, uL::AbstractVector, uR::AbstractVector)
    return rusanov_flux(eq, uL, uR)
end

