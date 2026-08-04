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
    residual!(du, state::SolutionState2D, eq, method)

Metric-aware 2D strong-form FR residual (Cartesian or curved) with capturing hooks.
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

    u_work = similar(state.u)
    preprocess_state!(u_work, method, state, eq)

    fill!(du, zero(T))

    # Face solution traces
    trW = zeros(T, Np, Neq, Nel)
    trE = zeros(T, Np, Neq, Nel)
    trS = zeros(T, Np, Neq, Nel)
    trN = zeros(T, Np, Neq, Nel)
    @inbounds for e in 1:Nel
        u_e = @view u_work[:, :, e, :]
        trW[:, :, e] = face_trace(u_e, ops, :west)
        trE[:, :, e] = face_trace(u_e, ops, :east)
        trS[:, :, e] = face_trace(u_e, ops, :south)
        trN[:, :, e] = face_trace(u_e, ops, :north)
    end

    # Continuous contravariant numerical fluxes on faces: F̃_hat, G̃_hat
    # Store as fhat_W/E (ξ-direction F̃) and fhat_S/N (η-direction G̃)
    fhat_W = zeros(T, Np, Neq, Nel)
    fhat_E = zeros(T, Np, Neq, Nel)
    fhat_S = zeros(T, Np, Neq, Nel)
    fhat_N = zeros(T, Np, Neq, Nel)

    flux_kind = state.scheme.flux

    # Interior vertical faces (shared F̃ = f̂_n * sJ with n from L→R = east of L)
    @inbounds for jy in 1:ny, jx in 1:(nx - 1)
        eL = element_index(mesh, jx, jy)
        eR = element_index(mesh, jx + 1, jy)
        for q in 1:Np
            u_m = collect(@view trE[q, :, eL])
            u_p = collect(@view trW[q, :, eR])
            # Use right-element-oriented continuous normal: east of L
            nnx, nny, sJ = met.nx_E[q, eL], met.ny_E[q, eL], met.sJ_E[q, eL]
            fh = interface_flux_n(eq, u_m, u_p, nnx, nny, flux_kind)
            for c in 1:Neq
                val = fh[c] * sJ
                fhat_E[q, c, eL] = val
                fhat_W[q, c, eR] = val
            end
        end
    end

    # Domain left/right
    @inbounds for jy in 1:ny
        eL = element_index(mesh, 1, jy)
        eR = element_index(mesh, nx, jy)
        for q in 1:Np
            if mesh.left_bc isa PeriodicBC
                u_m = collect(@view trE[q, :, eR])
                u_p = collect(@view trW[q, :, eL])
                nnx, nny, sJ = met.nx_E[q, eR], met.ny_E[q, eR], met.sJ_E[q, eR]
                fh = interface_flux_n(eq, u_m, u_p, nnx, nny, flux_kind)
                for c in 1:Neq
                    val = fh[c] * sJ
                    fhat_E[q, c, eR] = val
                    fhat_W[q, c, eL] = val
                end
            else
                # Transmissive BC: F̃ continuous uses +ξ geometric orientation
                # West face of eL: +ξ = -outward_W
                u_int = collect(@view trW[q, :, eL])
                sJ = met.sJ_W[q, eL]
                fh = interface_flux_n(
                    eq, u_int, u_int, -met.nx_W[q, eL], -met.ny_W[q, eL], flux_kind,
                )
                for c in 1:Neq
                    fhat_W[q, c, eL] = fh[c] * sJ
                end
                # East face of eR: +ξ = outward_E
                u_intR = collect(@view trE[q, :, eR])
                sJR = met.sJ_E[q, eR]
                fhR = interface_flux_n(
                    eq, u_intR, u_intR, met.nx_E[q, eR], met.ny_E[q, eR], flux_kind,
                )
                for c in 1:Neq
                    fhat_E[q, c, eR] = fhR[c] * sJR
                end
            end
        end
    end

    # Interior horizontal faces (shared G̃ = f̂_n * sJ, n = north of bottom)
    @inbounds for jy in 1:(ny - 1), jx in 1:nx
        eB = element_index(mesh, jx, jy)
        eT = element_index(mesh, jx, jy + 1)
        for q in 1:Np
            u_m = collect(@view trN[q, :, eB])
            u_p = collect(@view trS[q, :, eT])
            nnx, nny, sJ = met.nx_N[q, eB], met.ny_N[q, eB], met.sJ_N[q, eB]
            fh = interface_flux_n(eq, u_m, u_p, nnx, nny, flux_kind)
            for c in 1:Neq
                val = fh[c] * sJ
                fhat_N[q, c, eB] = val
                fhat_S[q, c, eT] = val
            end
        end
    end

    # Domain bottom/top
    @inbounds for jx in 1:nx
        eB = element_index(mesh, jx, 1)
        eT = element_index(mesh, jx, ny)
        for q in 1:Np
            if mesh.bottom_bc isa PeriodicBC
                u_m = collect(@view trN[q, :, eT])
                u_p = collect(@view trS[q, :, eB])
                nnx, nny, sJ = met.nx_N[q, eT], met.ny_N[q, eT], met.sJ_N[q, eT]
                fh = interface_flux_n(eq, u_m, u_p, nnx, nny, flux_kind)
                for c in 1:Neq
                    val = fh[c] * sJ
                    fhat_N[q, c, eT] = val
                    fhat_S[q, c, eB] = val
                end
            else
                # South of eB: +η = -outward_S
                u_int = collect(@view trS[q, :, eB])
                sJ = met.sJ_S[q, eB]
                fh = interface_flux_n(
                    eq, u_int, u_int, -met.nx_S[q, eB], -met.ny_S[q, eB], flux_kind,
                )
                for c in 1:Neq
                    fhat_S[q, c, eB] = fh[c] * sJ
                end
                # North of eT: +η = outward_N
                u_intT = collect(@view trN[q, :, eT])
                sJN = met.sJ_N[q, eT]
                fhT = interface_flux_n(
                    eq, u_intT, u_intT, met.nx_N[q, eT], met.ny_N[q, eT], flux_kind,
                )
                for c in 1:Neq
                    fhat_N[q, c, eT] = fhT[c] * sJN
                end
            end
        end
    end

    # Volume residual with metric fluxes
    @inbounds for e in 1:Nel
        # Build F̃, G̃ at all SPs
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

        # ξ-direction FR for each fixed j
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

        # η-direction FR for each fixed i
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
    apply_dissipation!(du, method, σ, u_work, state, eq)
    return du
end

function residual!(du, state::SolutionState2D, eq)
    return residual!(du, state, eq, NullCapturing())
end
