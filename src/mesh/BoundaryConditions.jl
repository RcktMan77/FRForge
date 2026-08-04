# Boundary condition types for 1D / 2D FR.

abstract type AbstractBC end

"""Periodic domain end (both ends must be PeriodicBC)."""
struct PeriodicBC <: AbstractBC end

"""Transmissive / zero-order extrapolation: ghost = interior trace."""
struct TransmissiveBC <: AbstractBC end

"""
Dirichlet: prescribed conserved state.

`u_func` may be:
- `t -> U` (1D style), or
- `(x, y, t) -> U` (2D face-local)
"""
struct DirichletBC{F} <: AbstractBC
    u_func::F
end

"""Slip / reflecting wall: reverse normal velocity, keep density and pressure."""
struct ReflectingBC <: AbstractBC end

"""
General exterior (ghost) state for mixed / problem-specific walls.

`ghost(u_int, nx_out, ny_out, x, y, t) -> AbstractVector` of conserved variables.
Outward unit normal is `(nx_out, ny_out)`.
"""
struct GhostStateBC{F} <: AbstractBC
    ghost::F
end

is_periodic(bc::AbstractBC) = bc isa PeriodicBC

"""
    reflect_conserved(U, nx, ny) -> U_ghost

Reflect normal velocity through unit outward normal `(nx, ny)`.
Kinetic energy (hence total energy) is unchanged.
"""
function reflect_conserved(U::AbstractVector{T}, nx, ny) where {T}
    ρ = max(U[1], eps(T))
    u = U[2] / ρ
    v = length(U) >= 3 ? U[3] / ρ : zero(T)
    un = u * T(nx) + v * T(ny)
    ug = u - 2 * un * T(nx)
    vg = v - 2 * un * T(ny)
    if length(U) == 2
        # 1D: U = (ρ, ρu) — rare; keep energy-less form
        return T[ρ, ρ * ug]
    elseif length(U) == 3
        # 1D Euler (ρ, ρu, E)
        return T[ρ, ρ * ug, U[3]]
    else
        # 2D Euler (ρ, ρu, ρv, E)
        return T[ρ, ρ * ug, ρ * vg, U[4]]
    end
end

"""
    exterior_state(bc, u_int, nx_out, ny_out, x, y, t) -> ghost conserved state

Resolve exterior state for any non-periodic BC. Used by 2D residual boundaries.
"""
function exterior_state(
    bc::TransmissiveBC,
    u_int::AbstractVector{T},
    nx_out,
    ny_out,
    x,
    y,
    t,
) where {T}
    return collect(u_int)
end

function exterior_state(
    bc::ReflectingBC,
    u_int::AbstractVector{T},
    nx_out,
    ny_out,
    x,
    y,
    t,
) where {T}
    return reflect_conserved(u_int, nx_out, ny_out)
end

function exterior_state(
    bc::DirichletBC,
    u_int::AbstractVector{T},
    nx_out,
    ny_out,
    x,
    y,
    t,
) where {T}
    f = bc.u_func
    ub = if applicable(f, x, y, t)
        f(x, y, t)
    else
        f(t)
    end
    Neq = length(u_int)
    out = Vector{T}(undef, Neq)
    @inbounds for c in 1:Neq
        out[c] = T(ub isa Number ? ub : ub[c])
    end
    return out
end

function exterior_state(
    bc::GhostStateBC,
    u_int::AbstractVector{T},
    nx_out,
    ny_out,
    x,
    y,
    t,
) where {T}
    ug = bc.ghost(u_int, nx_out, ny_out, x, y, t)
    Neq = length(u_int)
    out = Vector{T}(undef, Neq)
    @inbounds for c in 1:Neq
        out[c] = T(ug[c])
    end
    return out
end

function exterior_state(bc::AbstractBC, u_int, nx_out, ny_out, x, y, t)
    error("Unknown / unsupported BC type for exterior_state: $(typeof(bc))")
end
