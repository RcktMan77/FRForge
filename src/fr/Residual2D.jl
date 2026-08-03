# Tensor-product 2D FR residual on Cartesian quads.

"""
Extrapolate solution to one face of the reference element.

face ∈ (:west, :east, :south, :north) corresponding to ξ=-1, ξ=+1, η=-1, η=+1.
Returns matrix (Np, Neq) of face traces along the face SPs.
"""
function face_trace(u_e, ops::FROperators{T}, face::Symbol) where {T}
    # u_e[i,j,c]
    Np = size(u_e, 1)
    Neq = size(u_e, 3)
    out = zeros(T, Np, Neq)
    if face === :west  # ξ = -1, vary η = j
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
    elseif face === :south  # η = -1, vary ξ = i
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
    residual!(du, state::SolutionState2D, eq, method)

2D strong-form FR residual (tensor product) with staged capturing hooks:
preprocess → face traces / numerical flux → volume FR → sense + dissipate.

Hooks dispatch on `AbstractCapturingMethod` only (never concrete method names).
Face traces use Lagrange extrapolation (hybrid override path reserved for later).
"""
function residual!(
    du::AbstractArray{T,4},
    state::SolutionState2D{T,Neq},
    eq::AbstractEquation{Neq},
    method::AbstractCapturingMethod,
) where {T,Neq}
    mesh, ops = state.mesh, state.ops
    Np = n_points(ops)
    Nel = mesh.n_elements
    nx, ny = mesh.nx, mesh.ny

    u_work = similar(state.u)
    preprocess_state!(u_work, method, state, eq)

    fill!(du, zero(T))

    # --- Precompute all face traces per element ---
    # west/east: (Np, Neq, Nel), south/north: (Np, Neq, Nel)
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

    # Numerical fluxes on vertical faces: fhat_v[face_x, jy_face, jy_el_row_sp, c]
    # Vertical face between (jx,jy) and (jx+1,jy): normal +x
    # Store per-element left (west) and right (east) numerical fluxes: (Np, Neq, Nel)
    fhat_W = zeros(T, Np, Neq, Nel)
    fhat_E = zeros(T, Np, Neq, Nel)
    fhat_S = zeros(T, Np, Neq, Nel)
    fhat_N = zeros(T, Np, Neq, Nel)

    flux_kind = state.scheme.flux
    function flux_pair!(u_m, u_p, nx, ny)
        return interface_flux_n(eq, u_m, u_p, nx, ny, flux_kind)
    end

    # Interior vertical faces
    @inbounds for jy in 1:ny, jx in 1:(nx - 1)
        eL = element_index(mesh, jx, jy)
        eR = element_index(mesh, jx + 1, jy)
        for q in 1:Np  # along η
            u_m = @view trE[q, :, eL]
            u_p = @view trW[q, :, eR]
            fh = flux_pair!(collect(u_m), collect(u_p), 1, 0)
            for c in 1:Neq
                fhat_E[q, c, eL] = fh[c]
                fhat_W[q, c, eR] = fh[c]
            end
        end
    end

    # Domain left/right vertical faces
    @inbounds for jy in 1:ny
        eL = element_index(mesh, 1, jy)
        eR = element_index(mesh, nx, jy)
        for q in 1:Np
            if mesh.left_bc isa PeriodicBC
                u_m = collect(@view trE[q, :, eR])
                u_p = collect(@view trW[q, :, eL])
                fh = flux_pair!(u_m, u_p, 1, 0)
                for c in 1:Neq
                    fhat_W[q, c, eL] = fh[c]
                    fhat_E[q, c, eR] = fh[c]
                end
            else
                # Transmissive: ghost = interior
                u_int = collect(@view trW[q, :, eL])
                fh = flux_pair!(u_int, u_int, 1, 0)
                for c in 1:Neq
                    fhat_W[q, c, eL] = fh[c]
                end
                u_intR = collect(@view trE[q, :, eR])
                fhR = flux_pair!(u_intR, u_intR, 1, 0)
                for c in 1:Neq
                    fhat_E[q, c, eR] = fhR[c]
                end
            end
        end
    end

    # Interior horizontal faces
    @inbounds for jy in 1:(ny - 1), jx in 1:nx
        eB = element_index(mesh, jx, jy)
        eT = element_index(mesh, jx, jy + 1)
        for q in 1:Np  # along ξ
            u_m = collect(@view trN[q, :, eB])
            u_p = collect(@view trS[q, :, eT])
            fh = flux_pair!(u_m, u_p, 0, 1)
            for c in 1:Neq
                fhat_N[q, c, eB] = fh[c]
                fhat_S[q, c, eT] = fh[c]
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
                fh = flux_pair!(u_m, u_p, 0, 1)
                for c in 1:Neq
                    fhat_S[q, c, eB] = fh[c]
                    fhat_N[q, c, eT] = fh[c]
                end
            else
                u_int = collect(@view trS[q, :, eB])
                fh = flux_pair!(u_int, u_int, 0, 1)
                for c in 1:Neq
                    fhat_S[q, c, eB] = fh[c]
                end
                u_intT = collect(@view trN[q, :, eT])
                fhT = flux_pair!(u_intT, u_intT, 0, 1)
                for c in 1:Neq
                    fhat_N[q, c, eT] = fhT[c]
                end
            end
        end
    end

    # --- Volume residual ---
    D, gL, gR = ops.D, ops.gL_ξ, ops.gR_ξ
    @inbounds for e in 1:Nel
        jx, jy = element_coords(mesh, e)
        Jx, Jy = mesh.Jx[jx], mesh.Jy[jy]
        for j in 1:Np, i in 1:Np
            Uij = @view u_work[i, j, e, :]
            Fx = physical_flux_x(eq, Uij)
            Gy = physical_flux_y(eq, Uij)
            for c in 1:Neq
                # Need full lines for FR — compute per-line below
            end
        end
        # x-direction FR for each fixed j (η index)
        for j in 1:Np
            # F at SPs along i
            Fline = zeros(T, Np, Neq)
            for i in 1:Np
                Uij = @view u_work[i, j, e, :]
                fx = physical_flux_x(eq, Uij)
                for c in 1:Neq
                    Fline[i, c] = fx[c]
                end
            end
            for c in 1:Neq
                fL = zero(T)
                fR = zero(T)
                for i in 1:Np
                    fL += ops.ℓ_L[i] * Fline[i, c]
                    fR += ops.ℓ_R[i] * Fline[i, c]
                end
                for i in 1:Np
                    vol = zero(T)
                    for k in 1:Np
                        vol += D[i, k] * Fline[k, c]
                    end
                    corr =
                        (fhat_W[j, c, e] - fL) * gL[i] + (fhat_E[j, c, e] - fR) * gR[i]
                    du[i, j, e, c] -= (vol + corr) / Jx
                end
            end
        end
        # y-direction FR for each fixed i (ξ index)
        for i in 1:Np
            Gline = zeros(T, Np, Neq)
            for j in 1:Np
                Uij = @view u_work[i, j, e, :]
                gy = physical_flux_y(eq, Uij)
                for c in 1:Neq
                    Gline[j, c] = gy[c]
                end
            end
            for c in 1:Neq
                gS = zero(T)
                gN = zero(T)
                for j in 1:Np
                    gS += ops.ℓ_L[j] * Gline[j, c]
                    gN += ops.ℓ_R[j] * Gline[j, c]
                end
                for j in 1:Np
                    vol = zero(T)
                    for k in 1:Np
                        vol += D[j, k] * Gline[k, c]
                    end
                    corr =
                        (fhat_S[i, c, e] - gS) * gL[j] + (fhat_N[i, c, e] - gN) * gR[j]
                    du[i, j, e, c] -= (vol + corr) / Jy
                end
            end
        end
    end

    # Capturing: element sensor + dissipation (2D methods on Array{T,4})
    σ = zeros(T, Nel)
    sense!(σ, method, u_work, state, eq)
    apply_dissipation!(du, method, σ, u_work, state, eq)
    return du
end

function residual!(du, state::SolutionState2D, eq)
    return residual!(du, state, eq, NullCapturing())
end
