# 1D FR residual evaluation with BC ghosts and capturing hooks.

"""Interface numerical fluxes per element: fhat_L[e,c], fhat_R[e,c]."""
struct InterfaceFluxes{T}
    L::Matrix{T}  # (Nel, Neq)
    R::Matrix{T}
end

function _resolve_flux(eq, method, uL::AbstractVector, uR::AbstractVector)
    override = numerical_flux_method(method, eq, uL, uR)
    if override === nothing
        return numerical_flux(eq, uL, uR)
    end
    return override
end

"""
    compute_interface_fluxes(traces, mesh, eq, method, t) -> InterfaceFluxes

Apply BC ghosts at domain ends and compute one numerical flux per interface.
"""
function compute_interface_fluxes(
    traces::InterfaceTraces{T},
    mesh::Mesh1D{T},
    eq::AbstractEquation{Neq},
    method::AbstractCapturingMethod,
    t::T,
) where {T,Neq}
    Nel = mesh.n_elements
    fL = zeros(T, Nel, Neq)
    fR = zeros(T, Nel, Neq)

    # Helper to get flux from left/right state vectors
    function flux_at(u_minus::Vector{T}, u_plus::Vector{T})
        f = _resolve_flux(eq, method, u_minus, u_plus)
        return f isa AbstractVector ? collect(T, f) : T[T(f)]
    end

    # Interior interfaces between e and e+1
    for e in 1:(Nel - 1)
        u_m = Vector{T}(undef, Neq)
        u_p = Vector{T}(undef, Neq)
        for c in 1:Neq
            u_m[c] = traces.uR[e, c]
            u_p[c] = traces.uL[e + 1, c]
        end
        fhat = flux_at(u_m, u_p)
        for c in 1:Neq
            fR[e, c] = fhat[c]
            fL[e + 1, c] = fhat[c]
        end
    end

    # Domain left interface (element 1 left)
    u_m = Vector{T}(undef, Neq)
    u_p = Vector{T}(undef, Neq)
    for c in 1:Neq
        u_p[c] = traces.uL[1, c]
    end
    if mesh.left_bc isa PeriodicBC
        for c in 1:Neq
            u_m[c] = traces.uR[Nel, c]
        end
    elseif mesh.left_bc isa TransmissiveBC
        for c in 1:Neq
            u_m[c] = traces.uL[1, c]
        end
    elseif mesh.left_bc isa DirichletBC
        ub = mesh.left_bc.u_func(t)
        for c in 1:Neq
            u_m[c] = T(ub isa Number ? ub : ub[c])
        end
    else
        error("Unknown left BC type: $(typeof(mesh.left_bc))")
    end
    fhat_left = flux_at(u_m, u_p)
    for c in 1:Neq
        fL[1, c] = fhat_left[c]
        if mesh.left_bc isa PeriodicBC
            fR[Nel, c] = fhat_left[c]
        end
    end

    # Domain right interface (element Nel right) — skip if periodic (already set)
    if !(mesh.right_bc isa PeriodicBC)
        for c in 1:Neq
            u_m[c] = traces.uR[Nel, c]
        end
        if mesh.right_bc isa TransmissiveBC
            for c in 1:Neq
                u_p[c] = traces.uR[Nel, c]
            end
        elseif mesh.right_bc isa DirichletBC
            ub = mesh.right_bc.u_func(t)
            for c in 1:Neq
                u_p[c] = T(ub isa Number ? ub : ub[c])
            end
        else
            error("Unknown right BC type: $(typeof(mesh.right_bc))")
        end
        fhat_right = flux_at(u_m, u_p)
        for c in 1:Neq
            fR[Nel, c] = fhat_right[c]
        end
    end

    return InterfaceFluxes{T}(fL, fR)
end

"""Extrapolate discontinuous flux to left/right endpoints: sum_j f_j * ℓ_j(±1)."""
function extrapolate_flux(f::AbstractMatrix{T}, ℓ::AbstractVector{T}) where {T}
    Np, Neq = size(f)
    out = zeros(T, Neq)
    @inbounds for c in 1:Neq
        s = zero(T)
        for j in 1:Np
            s += ℓ[j] * f[j, c]
        end
        out[c] = s
    end
    return out
end

"""
    residual!(du, state, eq, method) -> du

Strong-form 1D FR residual with staged capturing hooks (Appendix A).
"""
function residual!(
    du::AbstractArray{T,3},
    state::SolutionState{T,Neq},
    eq::AbstractEquation{Neq},
    method::AbstractCapturingMethod,
) where {T,Neq}
    mesh, ops = state.mesh, state.ops
    Np, Nel, _ = size(state.u)

    u_work = similar(state.u)
    preprocess_state!(u_work, method, state, eq)

    traces = allocate_traces(Nel, Neq, T)
    extrapolate_interface!(traces, method, u_work, state, eq)

    fhat = compute_interface_fluxes(traces, mesh, eq, method, state.t)

    fill!(du, zero(T))
    @inbounds for e in 1:Nel
        U = @view u_work[:, e, :]
        f = physical_flux(eq, U)  # (Np, Neq)
        fL = extrapolate_flux(f, ops.ℓ_L)
        fR = extrapolate_flux(f, ops.ℓ_R)
        Je = mesh.J[e]
        for c in 1:Neq
            # vol = D * f[:, c]
            for j in 1:Np
                vol = zero(T)
                for k in 1:Np
                    vol += ops.D[j, k] * f[k, c]
                end
                corr =
                    (fhat.L[e, c] - fL[c]) * ops.gL_ξ[j] +
                    (fhat.R[e, c] - fR[c]) * ops.gR_ξ[j]
                du[j, e, c] = -(vol + corr) / Je
            end
        end
    end

    σ = zeros(T, Nel)
    sense!(σ, method, u_work, state, eq)
    apply_dissipation!(du, method, σ, u_work, state, eq)
    return du
end

residual!(du, state, eq) = residual!(du, state, eq, NullCapturing())

"""Check for non-finite entries in an array."""
function has_nonfinite(A)::Bool
    for x in A
        if !isfinite(x)
            return true
        end
    end
    return false
end
