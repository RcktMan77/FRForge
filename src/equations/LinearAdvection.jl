# 1D linear advection: u_t + a u_x = 0.

"""
    LinearAdvection1D{T}(a)

Constant-speed linear advection with wave speed `a`.
"""
struct LinearAdvection1D{T} <: AbstractEquation{1}
    a::T
end

LinearAdvection1D(a::Real) = LinearAdvection1D{typeof(float(a))}(float(a))

function physical_flux(eq::LinearAdvection1D{T}, u::AbstractVector) where {T}
    return T[eq.a * u[1]]
end

function physical_flux(eq::LinearAdvection1D{T}, u::Number) where {T}
    return eq.a * T(u)
end

"""Pointwise physical flux for a field of states (Np, Neq) → (Np, Neq)."""
function physical_flux!(out::AbstractVector{T}, eq::LinearAdvection1D{T}, u::AbstractVector) where {T}
    @inbounds out[1] = eq.a * u[1]
    return out
end

function physical_flux(eq::LinearAdvection1D{T}, U::AbstractMatrix) where {T}
    F = similar(U)
    physical_flux!(F, eq, U)
    return F
end

function max_wave_speed(eq::LinearAdvection1D, uL, uR)
    return abs(eq.a)
end

function max_wave_speed(eq::LinearAdvection1D, ::AbstractArray)
    return abs(eq.a)
end

"""Pure upwind numerical flux for linear advection."""
function numerical_flux(eq::LinearAdvection1D{T}, uL::AbstractVector, uR::AbstractVector) where {T}
    a = eq.a
    û = a >= 0 ? uL[1] : uR[1]
    return T[a * û]
end

function numerical_flux(eq::LinearAdvection1D{T}, uL::Number, uR::Number) where {T}
    a = eq.a
    û = a >= 0 ? T(uL) : T(uR)
    return a * û
end
