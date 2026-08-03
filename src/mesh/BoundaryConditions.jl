# Boundary condition types for 1D FR.

abstract type AbstractBC end

"""Periodic domain end (both ends must be PeriodicBC)."""
struct PeriodicBC <: AbstractBC end

"""Transmissive / zero-order extrapolation: ghost = interior trace."""
struct TransmissiveBC <: AbstractBC end

"""Dirichlet: prescribed conserved state as a function of time."""
struct DirichletBC{F} <: AbstractBC
    u_func::F  # u_func(t) -> AbstractVector of length Neq
end

is_periodic(bc::AbstractBC) = bc isa PeriodicBC
