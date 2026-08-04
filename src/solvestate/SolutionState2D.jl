# 2D solution state: u[ix, iy, e, c] tensor-product GL SPs.

"""
    SolutionState2D{T,Neq}

2D FR solution layout `u[i, j, e, c]` with tensor-product solution points.
Carries precomputed `metrics` for Cartesian or curved geometry.
"""
mutable struct SolutionState2D{T,Neq}
    u::Array{T,4}
    t::T
    p::Int
    mesh::Mesh2D{T}
    ops::FROperators{T}
    scheme::SchemeConfig
    metrics::MeshMetrics2D{T}
    residual_ws::Any  # ResidualWorkspace2D{T} | nothing (lazy)
end

function allocate_state(
    mesh::Mesh2D{T},
    ops::FROperators{T},
    ::Val{Neq};
    scheme::Union{Nothing,SchemeConfig} = nothing,
    metrics::Union{Nothing,MeshMetrics2D{T}} = nothing,
    wavy_amp::Union{Nothing,Real} = nothing,
) where {T,Neq}
    Np = n_points(ops)
    Nel = mesh.n_elements
    u = zeros(T, Np, Np, Nel, Neq)
    sch = something(
        scheme,
        SchemeConfig(; points = ops.points, flux = DEFAULT_SCHEME.flux, time = DEFAULT_SCHEME.time),
    )
    met = if metrics !== nothing
        metrics
    elseif wavy_amp !== nothing
        build_mesh_metrics_analytic_wavy(mesh, ops; amp = wavy_amp)
    else
        build_mesh_metrics(mesh, ops)
    end
    return SolutionState2D{T,Neq}(u, zero(T), ops.p, mesh, ops, sch, met, nothing)
end

function allocate_state(
    mesh::Mesh2D{T},
    ops::FROperators{T},
    neq::Int;
    scheme::Union{Nothing,SchemeConfig} = nothing,
    metrics::Union{Nothing,MeshMetrics2D{T}} = nothing,
    wavy_amp::Union{Nothing,Real} = nothing,
) where {T}
    return allocate_state(
        mesh, ops, Val(neq); scheme = scheme, metrics = metrics, wavy_amp = wavy_amp,
    )
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

"""Discrete integral of component `c` using |J| metric weights."""
function discrete_mass(state::SolutionState2D{T}, c::Int = 1) where {T}
    ops = state.ops
    met = state.metrics
    Np = n_points(ops)
    m = zero(T)
    for e in 1:state.mesh.n_elements
        for j in 1:Np, i in 1:Np
            m += ops.w[i] * ops.w[j] * abs(met.J[i, j, e]) * state.u[i, j, e, c]
        end
    end
    return m
end

function l2_error(state::SolutionState2D{T}, uexact, c::Int = 1) where {T}
    mesh, ops, met = state.mesh, state.ops, state.metrics
    Np = n_points(ops)
    acc = zero(T)
    for e in 1:mesh.n_elements
        for j in 1:Np, i in 1:Np
            x, y = physical_xy(mesh, e, ops.ξ[i], ops.ξ[j])
            val = uexact(x, y)
            exact = val isa Number ? T(val) : T(val[c])
            err = state.u[i, j, e, c] - exact
            acc += ops.w[i] * ops.w[j] * abs(met.J[i, j, e]) * err * err
        end
    end
    return sqrt(acc)
end
