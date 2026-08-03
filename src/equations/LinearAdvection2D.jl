# 2D linear advection: u_t + ax u_x + ay u_y = 0.

"""
    LinearAdvection2D{T}(ax, ay)

Constant-velocity 2D linear advection.
"""
struct LinearAdvection2D{T} <: AbstractEquation{1}
    ax::T
    ay::T
end

LinearAdvection2D(ax::Real, ay::Real) = LinearAdvection2D{typeof(float(ax))}(float(ax), float(ay))

"""x-direction physical flux F(u)."""
function physical_flux_x(eq::LinearAdvection2D{T}, u::AbstractVector) where {T}
    return T[eq.ax * u[1]]
end

function physical_flux_y(eq::LinearAdvection2D{T}, u::AbstractVector) where {T}
    return T[eq.ay * u[1]]
end

function max_wave_speed_n(eq::LinearAdvection2D, uL, uR, nx, ny)
    return abs(eq.ax * nx + eq.ay * ny)
end

function max_wave_speed(eq::LinearAdvection2D, ::AbstractArray)
    return hypot(eq.ax, eq.ay)
end

"""Upwind numerical flux in direction n=(nx,ny) for linear advection."""
function numerical_flux_n(eq::LinearAdvection2D{T}, uL::AbstractVector, uR::AbstractVector, nx, ny) where {T}
    an = eq.ax * T(nx) + eq.ay * T(ny)
    û = an >= 0 ? uL[1] : uR[1]
    return T[an * û]
end
