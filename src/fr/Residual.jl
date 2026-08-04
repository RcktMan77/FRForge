# 1D FR residual evaluation with BC ghosts and capturing hooks.

"""Interface numerical fluxes per element: fhat_L[e,c], fhat_R[e,c]."""
struct InterfaceFluxes{T}
    L::Matrix{T}  # (Nel, Neq)
    R::Matrix{T}
end

function _resolve_flux(eq, method, uL::AbstractVector, uR::AbstractVector, flux_kind::Symbol)
    override = numerical_flux_method(method, eq, uL, uR)
    if override === nothing
        return interface_flux(eq, uL, uR, flux_kind)
    end
    return override
end

"""
    compute_interface_fluxes!(ws, traces, mesh, eq, method, t; flux_kind=:rusanov) -> InterfaceFluxes

Apply BC ghosts at domain ends and compute one numerical flux per interface.
Writes into workspace face buffers (`ws.fL` / `ws.fR`). Prefer this over the
allocating `compute_interface_fluxes` wrapper.
"""
function compute_interface_fluxes!(
    ws::ResidualWorkspace1D{T},
    traces::InterfaceTraces{T},
    mesh::Mesh1D{T},
    eq::AbstractEquation{Neq},
    method::AbstractCapturingMethod,
    t::T;
    flux_kind::Symbol = :rusanov,
) where {T,Neq}
    Nel = mesh.n_elements
    fL, fR = ws.fL, ws.fR
    u_m, u_p = ws.u_m, ws.u_p

    function flux_at!(dest::AbstractVector{T}, u_minus::Vector{T}, u_plus::Vector{T})
        f = _resolve_flux(eq, method, u_minus, u_plus, flux_kind)
        if f isa AbstractVector
            @inbounds for c in 1:Neq
                dest[c] = T(f[c])
            end
        else
            dest[1] = T(f)
        end
        return dest
    end
    ftmp = Vector{T}(undef, Neq)

    # Interior interfaces between e and e+1
    for e in 1:(Nel - 1)
        for c in 1:Neq
            u_m[c] = traces.uR[e, c]
            u_p[c] = traces.uL[e + 1, c]
        end
        flux_at!(ftmp, u_m, u_p)
        for c in 1:Neq
            fR[e, c] = ftmp[c]
            fL[e + 1, c] = ftmp[c]
        end
    end

    # Domain left interface (element 1 left)
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
    flux_at!(ftmp, u_m, u_p)
    for c in 1:Neq
        fL[1, c] = ftmp[c]
        if mesh.left_bc isa PeriodicBC
            fR[Nel, c] = ftmp[c]
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
        flux_at!(ftmp, u_m, u_p)
        for c in 1:Neq
            fR[Nel, c] = ftmp[c]
        end
    end

    return InterfaceFluxes{T}(fL, fR)
end

"""
    compute_interface_fluxes(traces, mesh, eq, method, t; flux_kind=:rusanov) -> InterfaceFluxes

Allocating convenience wrapper (creates a temporary workspace). Prefer
`compute_interface_fluxes!` from the residual path.
"""
function compute_interface_fluxes(
    traces::InterfaceTraces{T},
    mesh::Mesh1D{T},
    eq::AbstractEquation{Neq},
    method::AbstractCapturingMethod,
    t::T;
    flux_kind::Symbol = :rusanov,
) where {T,Neq}
    Nel = mesh.n_elements
    ws = ResidualWorkspace1D(T, 1, Nel, Neq)  # Np unused for fluxes
    return compute_interface_fluxes!(ws, traces, mesh, eq, method, t; flux_kind = flux_kind)
end

"""Extrapolate discontinuous flux to left/right endpoints: sum_j f_j * ℓ_j(±1)."""
function extrapolate_flux!(out::AbstractVector{T}, f::AbstractMatrix{T}, ℓ::AbstractVector{T}) where {T}
    Np, Neq = size(f)
    length(out) == Neq || throw(DimensionMismatch("out length $(length(out)) != Neq $Neq"))
    @inbounds for c in 1:Neq
        s = zero(T)
        for j in 1:Np
            s += ℓ[j] * f[j, c]
        end
        out[c] = s
    end
    return out
end

function extrapolate_flux(f::AbstractMatrix{T}, ℓ::AbstractVector{T}) where {T}
    out = Vector{T}(undef, size(f, 2))
    return extrapolate_flux!(out, f, ℓ)
end

"""
    residual!(du, state, eq, method) -> du

Strong-form 1D FR residual with staged capturing hooks (see design residual pipeline).
Reuses `state` residual workspace (`u_work`, traces, interface fluxes, `σ`).
"""
function residual!(
    du::AbstractArray{T,3},
    state::SolutionState{T,Neq},
    eq::AbstractEquation{Neq},
    method::AbstractCapturingMethod,
) where {T,Neq}
    mesh, ops = state.mesh, state.ops
    Np, Nel, _ = size(state.u)

    ws = ensure_residual_workspace!(state)
    u_work = ws.u_work
    preprocess_state!(u_work, method, state, eq)

    traces = ws.traces
    extrapolate_interface!(traces, method, u_work, state, eq)

    fhat = compute_interface_fluxes!(
        ws,
        traces,
        mesh,
        eq,
        method,
        state.t;
        flux_kind = state.scheme.flux,
    )

    fill!(du, zero(T))
    f = ws.f_vol
    fL, fR = ws.f_end_L, ws.f_end_R
    @inbounds for e in 1:Nel
        U = @view u_work[:, e, :]
        physical_flux!(f, eq, U)  # (Np, Neq) into workspace
        extrapolate_flux!(fL, f, ops.ℓ_L)
        extrapolate_flux!(fR, f, ops.ℓ_R)
        Je = mesh.J[e]
        for c in 1:Neq
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

    σ = ws.σ
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
