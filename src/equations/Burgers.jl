# 1D inviscid Burgers: u_t + (u²/2)_x = 0.

"""
    Burgers1D

Inviscid Burgers equation (scalar nonlinear conservation law).
Physical flux F(u) = u²/2.
"""
struct Burgers1D <: AbstractEquation{1} end

function physical_flux(::Burgers1D, u::AbstractVector{T}) where {T}
    return T[T(0.5) * u[1] * u[1]]
end

function physical_flux(::Burgers1D, u::Number)
    return 0.5 * float(u)^2
end

function physical_flux!(out::AbstractVector{T}, ::Burgers1D, u::AbstractVector) where {T}
    @inbounds out[1] = T(0.5) * u[1] * u[1]
    return out
end

function physical_flux(::Burgers1D, U::AbstractMatrix{T}) where {T}
    F = similar(U)
    physical_flux!(F, Burgers1D(), U)
    return F
end

function max_wave_speed(::Burgers1D, uL::AbstractVector, uR::AbstractVector)
    return max(abs(uL[1]), abs(uR[1]))
end

function max_wave_speed(::Burgers1D, uL::Number, uR::Number)
    return max(abs(uL), abs(uR))
end

"""Global max |u| over the solution array (for CFL)."""
function max_wave_speed(::Burgers1D, U::AbstractArray{T}) where {T}
    m = zero(T)
    @inbounds for x in U
        m = max(m, abs(x))
    end
    return max(m, eps(T))
end

"""Rusanov (local Lax–Friedrichs) numerical flux for Burgers."""
function numerical_flux(eq::Burgers1D, uL::AbstractVector{T}, uR::AbstractVector{T}) where {T}
    return rusanov_flux(eq, uL, uR)
end
