# 2D solution state: u[ix, iy, e, c] tensor-product GL SPs.

"""
    SolutionState2D{T,Neq}

2D FR solution layout `u[i, j, e, c]` with tensor-product solution points.
"""
mutable struct SolutionState2D{T,Neq}
    u::Array{T,4}
    t::T
    p::Int
    mesh::Mesh2D{T}
    ops::FROperators{T}
    scheme::SchemeConfig
end

function allocate_state(
    mesh::Mesh2D{T},
    ops::FROperators{T},
    ::Val{Neq};
    scheme::Union{Nothing,SchemeConfig}=nothing,
) where {T,Neq}
    Np = n_points(ops)
    Nel = mesh.n_elements
    u = zeros(T, Np, Np, Nel, Neq)
    sch = something(
        scheme,
        SchemeConfig(; points=ops.points, flux=DEFAULT_SCHEME.flux, time=DEFAULT_SCHEME.time),
    )
    return SolutionState2D{T,Neq}(u, zero(T), ops.p, mesh, ops, sch)
end

function allocate_state(
    mesh::Mesh2D{T},
    ops::FROperators{T},
    neq::Int;
    scheme::Union{Nothing,SchemeConfig}=nothing,
) where {T}
    return allocate_state(mesh, ops, Val(neq); scheme=scheme)
end

"""
    set_initial_condition!(state::SolutionState2D, f)

`f(x, y)` returns a Number (Neq=1) or AbstractVector of length Neq.
"""
function set_initial_condition!(state::SolutionState2D{T,Neq}, f) where {T,Neq}
    mesh, ops = state.mesh, state.ops
    Np = n_points(ops)
    Nel = mesh.n_elements
    for e in 1:Nel
        for jy in 1:Np, jx in 1:Np
            x, y = physical_xy(mesh, e, ops.ξ[jx], ops.ξ[jy])
            val = f(x, y)
            if Neq == 1
                state.u[jx, jy, e, 1] = T(val isa Number ? val : val[1])
            else
                for c in 1:Neq
                    state.u[jx, jy, e, c] = T(val[c])
                end
            end
        end
    end
    state.t = zero(T)
    return state
end

"""Discrete integral of component `c` (GL tensor weights × |J| = Jx*Jy * 4 area factor).

Reference weights w on [-1,1]; physical ∫ = Σ w_i w_j Jx Jy u.
"""
function discrete_mass(state::SolutionState2D{T}, c::Int=1) where {T}
    mesh, ops = state.mesh, state.ops
    Np = n_points(ops)
    m = zero(T)
    for e in 1:mesh.n_elements
        jx, jy = element_coords(mesh, e)
        Jx, Jy = mesh.Jx[jx], mesh.Jy[jy]
        for j in 1:Np, i in 1:Np
            m += ops.w[i] * ops.w[j] * Jx * Jy * state.u[i, j, e, c]
        end
    end
    return m
end

function l2_error(state::SolutionState2D{T}, uexact, c::Int=1) where {T}
    mesh, ops = state.mesh, state.ops
    Np = n_points(ops)
    acc = zero(T)
    for e in 1:mesh.n_elements
        jx, jy = element_coords(mesh, e)
        Jx, Jy = mesh.Jx[jx], mesh.Jy[jy]
        for j in 1:Np, i in 1:Np
            x, y = physical_xy(mesh, e, ops.ξ[i], ops.ξ[j])
            val = uexact(x, y)
            exact = val isa Number ? T(val) : T(val[c])
            err = state.u[i, j, e, c] - exact
            acc += ops.w[i] * ops.w[j] * Jx * Jy * err * err
        end
    end
    return sqrt(acc)
end

