# Abstract equation interface.

"""
    AbstractEquation{Neq}

Hyperbolic conservation law with `Neq` conserved components.
"""
abstract type AbstractEquation{Neq} end

n_equations(::AbstractEquation{Neq}) where {Neq} = Neq

"""Physical flux F(u) at a single state; returns vector of length Neq."""
function physical_flux end

"""Maximum wave speed for CFL (scalar >= 0)."""
function max_wave_speed end

"""Default numerical flux at an interface from left/right states."""
function numerical_flux end
