# Tensor-product 2D FR residual with metric terms (Cartesian or curved quads).

"""
Extrapolate solution to one face of the reference element.

face ∈ (:west, :east, :south, :north) corresponding to ξ=-1, ξ=+1, η=-1, η=+1.
Returns matrix (Np, Neq) of face traces along the face SPs.
"""
function face_trace(u_e, ops::FROperators{T}, face::Symbol) where {T}
    Np = size(u_e, 1)
    Neq = size(u_e, 3)
    out = zeros(T, Np, Neq)
    if face === :west
        @inbounds for j in 1:Np, c in 1:Neq
            s = zero(T)
            for i in 1:Np
                s += ops.ℓ_L[i] * u_e[i, j, c]
            end
            out[j, c] = s
        end
    elseif face === :east
        @inbounds for j in 1:Np, c in 1:Neq
            s = zero(T)
            for i in 1:Np
                s += ops.ℓ_R[i] * u_e[i, j, c]
            end
            out[j, c] = s
        end
    elseif face === :south
        @inbounds for i in 1:Np, c in 1:Neq
            s = zero(T)
            for j in 1:Np
                s += ops.ℓ_L[j] * u_e[i, j, c]
            end
            out[i, c] = s
        end
    elseif face === :north
        @inbounds for i in 1:Np, c in 1:Neq
            s = zero(T)
            for j in 1:Np
                s += ops.ℓ_R[j] * u_e[i, j, c]
            end
            out[i, c] = s
        end
    else
        error("unknown face $face")
    end
    return out
end

"""
Contravariant fluxes at a solution point:
  F̃ =  y_η F_x - x_η F_y
  G̃ = -y_ξ F_x + x_ξ F_y
"""
@inline function contravariant_fluxes(
    Fx::AbstractVector{T},
    Gy::AbstractVector{T},
    x_ξ::T,
    x_η::T,
    y_ξ::T,
    y_η::T,
) where {T}
    Neq = length(Fx)
    Ft = Vector{T}(undef, Neq)
    Gt = Vector{T}(undef, Neq)
    @inbounds for c in 1:Neq
        Ft[c] = y_η * Fx[c] - x_η * Gy[c]
        Gt[c] = -y_ξ * Fx[c] + x_ξ * Gy[c]
    end
    return Ft, Gt
end

"""
Fill continuous contravariant face flux (F̃ or G̃) from interior state + BC.

`into_domain`: true when the reference +direction points into the domain
(west/south faces). Then L=ghost, R=interior along +ref.
"""
function boundary_fhat!(
    fhat_face::AbstractArray{T,3},
    q::Int,
    e::Int,
    eq,
    bc::AbstractBC,
    u_int::AbstractVector{T},
    nx_out::T,
    ny_out::T,
    sJ::T,
    x::T,
    y::T,
    t::T,
    flux_kind::Symbol,
    into_domain::Bool,
) where {T}
    ug = exterior_state(bc, u_int, nx_out, ny_out, x, y, t)
    if into_domain
        # +ref = -outward: L=ghost, R=int, n_+ref = -n_out
        fh = interface_flux_n(eq, ug, u_int, -nx_out, -ny_out, flux_kind)
    else
        # +ref = +outward: L=int, R=ghost, n_+ref = n_out
        fh = interface_flux_n(eq, u_int, ug, nx_out, ny_out, flux_kind)
    end
    @inbounds for c in 1:length(fh)
        fhat_face[q, c, e] = fh[c] * sJ
    end
    return nothing
end

"""
    residual!(du, state::SolutionState2D, eq, method)

Metric-aware 2D strong-form FR residual (Cartesian or curved) with capturing hooks.
Supports Periodic / Transmissive / Reflecting / Dirichlet / GhostState BCs and
optional solid-element masks (forward-facing step etc.).
"""
function residual!(
    du::AbstractArray{T,4},
    state::SolutionState2D{T,Neq},
    eq::AbstractEquation{Neq},
    method::AbstractCapturingMethod,
) where {T,Neq}
    mesh, ops, met = state.mesh, state.ops, state.metrics
    Np = n_points(ops)
    Nel = mesh.n_elements
    nx, ny = mesh.nx, mesh.ny
    D, gL, gR = ops.D, ops.gL_ξ, ops.gR_ξ
    ℓ_L, ℓ_R = ops.ℓ_L, ops.ℓ_R
    t = state.t
    has_solid = mesh.solid !== nothing

    u_work = similar(state.u)
    preprocess_state!(u_work, method, state, eq)

    fill!(du, zero(T))

    # Face solution traces
    trW = zeros(T, Np, Neq, Nel)
    trE = zeros(T, Np, Neq, Nel)
    trS = zeros(T, Np, Neq, Nel)
    trN = zeros(T, Np, Neq, Nel)
    @inbounds for e in 1:Nel
        has_solid && is_solid(mesh, e) && continue
        u_e = @view u_work[:, :, e, :]
        trW[:, :, e] = face_trace(u_e, ops, :west)
        trE[:, :, e] = face_trace(u_e, ops, :east)
        trS[:, :, e] = face_trace(u_e, ops, :south)
        trN[:, :, e] = face_trace(u_e, ops, :north)
    end

    # Continuous contravariant numerical fluxes on faces: F̃_hat, G̃_hat
    fhat_W = zeros(T, Np, Neq, Nel)
    fhat_E = zeros(T, Np, Neq, Nel)
    fhat_S = zeros(T, Np, Neq, Nel)
    fhat_N = zeros(T, Np, Neq, Nel)

    flux_kind = state.scheme.flux
    wall = ReflectingBC()

    # Interior vertical faces
    @inbounds for jy in 1:ny, jx in 1:(nx - 1)
        eL = element_index(mesh, jx, jy)
        eR = element_index(mesh, jx + 1, jy)
        sL = has_solid && is_solid(mesh, eL)
        sR = has_solid && is_solid(mesh, eR)
        (sL && sR) && continue
        for q in 1:Np
            nnx, nny, sJ = met.nx_E[q, eL], met.ny_E[q, eL], met.sJ_E[q, eL]
            if sL && !sR
                # Solid on left → wall on west of fluid eR
                u_int = collect(@view trW[q, :, eR])
                x, y = physical_xy(mesh, eR, -one(T), ops.ξ[q])
                boundary_fhat!(
                    fhat_W, q, eR, eq, wall, u_int,
                    met.nx_W[q, eR], met.ny_W[q, eR], met.sJ_W[q, eR],
                    x, y, t, flux_kind, true,
                )
            elseif sR && !sL
                # Solid on right → wall on east of fluid eL
                u_int = collect(@view trE[q, :, eL])
                x, y = physical_xy(mesh, eL, one(T), ops.ξ[q])
                boundary_fhat!(
                    fhat_E, q, eL, eq, wall, u_int, nnx, nny, sJ,
                    x, y, t, flux_kind, false,
                )
            else
                u_m = collect(@view trE[q, :, eL])
                u_p = collect(@view trW[q, :, eR])
                fh = interface_flux_n(eq, u_m, u_p, nnx, nny, flux_kind)
                for c in 1:Neq
                    val = fh[c] * sJ
                    fhat_E[q, c, eL] = val
                    fhat_W[q, c, eR] = val
                end
            end
        end
    end

    # Domain left/right
    @inbounds for jy in 1:ny
        eL = element_index(mesh, 1, jy)
        eR = element_index(mesh, nx, jy)
        for q in 1:Np
            if mesh.left_bc isa PeriodicBC
                if !(has_solid && (is_solid(mesh, eL) || is_solid(mesh, eR)))
                    u_m = collect(@view trE[q, :, eR])
                    u_p = collect(@view trW[q, :, eL])
                    nnx, nny, sJ = met.nx_E[q, eR], met.ny_E[q, eR], met.sJ_E[q, eR]
                    fh = interface_flux_n(eq, u_m, u_p, nnx, nny, flux_kind)
                    for c in 1:Neq
                        val = fh[c] * sJ
                        fhat_E[q, c, eR] = val
                        fhat_W[q, c, eL] = val
                    end
                end
            else
                if !(has_solid && is_solid(mesh, eL))
                    u_int = collect(@view trW[q, :, eL])
                    x, y = physical_xy(mesh, eL, -one(T), ops.ξ[q])
                    boundary_fhat!(
                        fhat_W, q, eL, eq, mesh.left_bc, u_int,
                        met.nx_W[q, eL], met.ny_W[q, eL], met.sJ_W[q, eL],
                        x, y, t, flux_kind, true,
                    )
                end
                if !(has_solid && is_solid(mesh, eR))
                    u_intR = collect(@view trE[q, :, eR])
                    xR, yR = physical_xy(mesh, eR, one(T), ops.ξ[q])
                    boundary_fhat!(
                        fhat_E, q, eR, eq, mesh.right_bc, u_intR,
                        met.nx_E[q, eR], met.ny_E[q, eR], met.sJ_E[q, eR],
                        xR, yR, t, flux_kind, false,
                    )
                end
            end
        end
    end

    # Interior horizontal faces
    @inbounds for jy in 1:(ny - 1), jx in 1:nx
        eB = element_index(mesh, jx, jy)
        eT = element_index(mesh, jx, jy + 1)
        sB = has_solid && is_solid(mesh, eB)
        sT = has_solid && is_solid(mesh, eT)
        (sB && sT) && continue
        for q in 1:Np
            nnx, nny, sJ = met.nx_N[q, eB], met.ny_N[q, eB], met.sJ_N[q, eB]
            if sB && !sT
                u_int = collect(@view trS[q, :, eT])
                x, y = physical_xy(mesh, eT, ops.ξ[q], -one(T))
                boundary_fhat!(
                    fhat_S, q, eT, eq, wall, u_int,
                    met.nx_S[q, eT], met.ny_S[q, eT], met.sJ_S[q, eT],
                    x, y, t, flux_kind, true,
                )
            elseif sT && !sB
                u_int = collect(@view trN[q, :, eB])
                x, y = physical_xy(mesh, eB, ops.ξ[q], one(T))
                boundary_fhat!(
                    fhat_N, q, eB, eq, wall, u_int, nnx, nny, sJ,
                    x, y, t, flux_kind, false,
                )
            else
                u_m = collect(@view trN[q, :, eB])
                u_p = collect(@view trS[q, :, eT])
                fh = interface_flux_n(eq, u_m, u_p, nnx, nny, flux_kind)
                for c in 1:Neq
                    val = fh[c] * sJ
                    fhat_N[q, c, eB] = val
                    fhat_S[q, c, eT] = val
                end
            end
        end
    end

    # Domain bottom/top
    @inbounds for jx in 1:nx
        eB = element_index(mesh, jx, 1)
        eT = element_index(mesh, jx, ny)
        for q in 1:Np
            if mesh.bottom_bc isa PeriodicBC
                if !(has_solid && (is_solid(mesh, eB) || is_solid(mesh, eT)))
                    u_m = collect(@view trN[q, :, eT])
                    u_p = collect(@view trS[q, :, eB])
                    nnx, nny, sJ = met.nx_N[q, eT], met.ny_N[q, eT], met.sJ_N[q, eT]
                    fh = interface_flux_n(eq, u_m, u_p, nnx, nny, flux_kind)
                    for c in 1:Neq
                        val = fh[c] * sJ
                        fhat_N[q, c, eT] = val
                        fhat_S[q, c, eB] = val
                    end
                end
            else
                if !(has_solid && is_solid(mesh, eB))
                    u_int = collect(@view trS[q, :, eB])
                    x, y = physical_xy(mesh, eB, ops.ξ[q], -one(T))
                    boundary_fhat!(
                        fhat_S, q, eB, eq, mesh.bottom_bc, u_int,
                        met.nx_S[q, eB], met.ny_S[q, eB], met.sJ_S[q, eB],
                        x, y, t, flux_kind, true,
                    )
                end
                if !(has_solid && is_solid(mesh, eT))
                    u_intT = collect(@view trN[q, :, eT])
                    xT, yT = physical_xy(mesh, eT, ops.ξ[q], one(T))
                    boundary_fhat!(
                        fhat_N, q, eT, eq, mesh.top_bc, u_intT,
                        met.nx_N[q, eT], met.ny_N[q, eT], met.sJ_N[q, eT],
                        xT, yT, t, flux_kind, false,
                    )
                end
            end
        end
    end

    # Volume residual with metric fluxes (skip solid elements)
    @inbounds for e in 1:Nel
        has_solid && is_solid(mesh, e) && continue
        Ft = zeros(T, Np, Np, Neq)
        Gt = zeros(T, Np, Np, Neq)
        for j in 1:Np, i in 1:Np
            Uij = @view u_work[i, j, e, :]
            Fx = physical_flux_x(eq, Uij)
            Gy = physical_flux_y(eq, Uij)
            xξ, xη = met.x_ξ[i, j, e], met.x_η[i, j, e]
            yξ, yη = met.y_ξ[i, j, e], met.y_η[i, j, e]
            for c in 1:Neq
                Ft[i, j, c] = yη * Fx[c] - xη * Gy[c]
                Gt[i, j, c] = -yξ * Fx[c] + xξ * Gy[c]
            end
        end

        for j in 1:Np
            for c in 1:Neq
                fL = zero(T)
                fR = zero(T)
                for i in 1:Np
                    fL += ℓ_L[i] * Ft[i, j, c]
                    fR += ℓ_R[i] * Ft[i, j, c]
                end
                for i in 1:Np
                    vol = zero(T)
                    for k in 1:Np
                        vol += D[i, k] * Ft[k, j, c]
                    end
                    corr =
                        (fhat_W[j, c, e] - fL) * gL[i] + (fhat_E[j, c, e] - fR) * gR[i]
                    Je = met.J[i, j, e]
                    du[i, j, e, c] -= (vol + corr) / Je
                end
            end
        end

        for i in 1:Np
            for c in 1:Neq
                gS = zero(T)
                gN = zero(T)
                for j in 1:Np
                    gS += ℓ_L[j] * Gt[i, j, c]
                    gN += ℓ_R[j] * Gt[i, j, c]
                end
                for j in 1:Np
                    vol = zero(T)
                    for k in 1:Np
                        vol += D[j, k] * Gt[i, k, c]
                    end
                    corr =
                        (fhat_S[i, c, e] - gS) * gL[j] + (fhat_N[i, c, e] - gN) * gR[j]
                    Je = met.J[i, j, e]
                    du[i, j, e, c] -= (vol + corr) / Je
                end
            end
        end
    end

    σ = zeros(T, Nel)
    sense!(σ, method, u_work, state, eq)
    if has_solid
        @inbounds for e in 1:Nel
            is_solid(mesh, e) && (σ[e] = zero(T))
        end
    end
    apply_dissipation!(du, method, σ, u_work, state, eq)
    if has_solid
        @inbounds for e in 1:Nel
            if is_solid(mesh, e)
                du[:, :, e, :] .= zero(T)
            end
        end
    end
    return du
end

function residual!(du, state::SolutionState2D, eq)
    return residual!(du, state, eq, NullCapturing())
end
