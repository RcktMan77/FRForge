# 1D solution state: u[Np, Nel, Neq]

"""
    SolutionState{T,Neq}

1D FR solution: conserved variables at solution points, time, mesh, operators.
Layout: `u[j, e, c]` with j solution point, e element, c equation component.
"""
mutable struct SolutionState{T,Neq}
    u::Array{T,3}
    t::T
    p::Int
    mesh::Mesh1D{T}
    ops::FROperators{T}
    scheme::SchemeConfig
end

"""
    allocate_state(mesh, ops, ::Val{Neq}; scheme) -> SolutionState

Allocate zero-initialized state for `Neq` conserved components.
`scheme` defaults to DEFAULT_SCHEME with `points` taken from `ops.points`.
"""
function allocate_state(
    mesh::Mesh1D{T},
    ops::FROperators{T},
    ::Val{Neq};
    scheme::Union{Nothing,SchemeConfig}=nothing,
) where {T,Neq}
    Np = n_points(ops)
    Nel = mesh.n_elements
    u = zeros(T, Np, Nel, Neq)
    sch = something(
        scheme,
        SchemeConfig(; points=ops.points, flux=DEFAULT_SCHEME.flux, time=DEFAULT_SCHEME.time),
    )
    return SolutionState{T,Neq}(u, zero(T), ops.p, mesh, ops, sch)
end

function allocate_state(
    mesh::Mesh1D{T},
    ops::FROperators{T},
    neq::Int;
    scheme::Union{Nothing,SchemeConfig}=nothing,
) where {T}
    return allocate_state(mesh, ops, Val(neq); scheme=scheme)
end

"""
    set_initial_condition!(state, u0!)

`u0!(out, x, e)` fills `out` (length Neq) at physical coordinate `x` in element `e`,
or for scalar equations `u0(x)` may be used via a wrapper.

Simpler API: `set_initial_condition!(state, f)` where `f(x) -> Number` or
`f(x) -> AbstractVector` of length Neq.
"""
function set_initial_condition!(state::SolutionState{T,Neq}, f) where {T,Neq}
    mesh, ops = state.mesh, state.ops
    Np, Nel = size(state.u, 1), size(state.u, 2)
    for e in 1:Nel
        xs = physical_coords(mesh, ops, e)
        for j in 1:Np
            val = f(xs[j])
            if Neq == 1
                state.u[j, e, 1] = T(val isa Number ? val : val[1])
            else
                for c in 1:Neq
                    state.u[j, e, c] = T(val[c])
                end
            end
        end
    end
    state.t = zero(T)
    return state
end

"""Discrete mass (integral) of component `c` using GL quadrature."""
function discrete_mass(state::SolutionState{T}, c::Int=1) where {T}
    mesh, ops = state.mesh, state.ops
    Np, Nel = size(state.u, 1), size(state.u, 2)
    m = zero(T)
    for e in 1:Nel
        for j in 1:Np
            m += ops.w[j] * mesh.J[e] * state.u[j, e, c]
        end
    end
    return m
end

"""L2 norm of error vs exact function `uexact(x)` for component `c` (scalar exact)."""
function l2_error(state::SolutionState{T}, uexact, c::Int=1) where {T}
    mesh, ops = state.mesh, state.ops
    Np, Nel = size(state.u, 1), size(state.u, 2)
    acc = zero(T)
    for e in 1:Nel
        xs = physical_coords(mesh, ops, e)
        for j in 1:Np
            val = uexact(xs[j])
            exact = val isa Number ? T(val) : T(val[c])
            err = state.u[j, e, c] - exact
            acc += ops.w[j] * mesh.J[e] * err * err
        end
    end
    return sqrt(acc)
end

"""Combined L2 error over all conserved components: √(Σ_c ‖U_c−Uex_c‖²)."""
function l2_error_all(state::SolutionState{T,Neq}, uexact) where {T,Neq}
    s = zero(T)
    for c in 1:Neq
        e = l2_error(state, uexact, c)
        s += e * e
    end
    return sqrt(s)
end