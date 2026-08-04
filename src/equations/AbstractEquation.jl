# Abstract equation interface.

"""
    AbstractEquation{Neq}

Hyperbolic conservation law with `Neq` conserved components.
"""
abstract type AbstractEquation{Neq} end

n_equations(::AbstractEquation{Neq}) where {Neq} = Neq

"""Physical flux F(u) at a single state; returns vector of length Neq."""
function physical_flux end

"""In-place physical flux for a single state vector (default: allocate then copy)."""
function physical_flux!(out::AbstractVector, eq::AbstractEquation, U::AbstractVector)
    f = physical_flux(eq, U)
    @inbounds for c in 1:length(out)
        out[c] = f[c]
    end
    return out
end

"""In-place physical flux for SP field `U` (Np, Neq) into `F` (Np, Neq)."""
function physical_flux!(F::AbstractMatrix{T}, eq::AbstractEquation, U::AbstractMatrix{T}) where {T}
    Np, Neq = size(U)
    size(F) == (Np, Neq) || throw(DimensionMismatch("F size $(size(F)) != U size $(size(U))"))
    @inbounds for j in 1:Np
        physical_flux!(view(F, j, :), eq, view(U, j, :))
    end
    return F
end

# Optional in-place 2D physical fluxes (default: allocate then copy)
function physical_flux_x!(out::AbstractVector, eq::AbstractEquation, U)
    Fx = physical_flux_x(eq, U)
    @inbounds for c in 1:length(out)
        out[c] = Fx[c]
    end
    return out
end

function physical_flux_y!(out::AbstractVector, eq::AbstractEquation, U)
    Gy = physical_flux_y(eq, U)
    @inbounds for c in 1:length(out)
        out[c] = Gy[c]
    end
    return out
end

function numerical_flux_n!(out::AbstractVector, eq::AbstractEquation, uL, uR, nx, ny)
    fh = numerical_flux_n(eq, uL, uR, nx, ny)
    @inbounds for c in 1:length(out)
        out[c] = fh[c]
    end
    return out
end

"""Maximum wave speed for CFL (scalar >= 0)."""
function max_wave_speed end

"""Default numerical flux at an interface from left/right states."""
function numerical_flux end
