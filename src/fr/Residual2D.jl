# Tensor-product 2D FR residual with metric terms (Cartesian or curved quads).
# Phase 4: buffer reuse / in-place fluxes — numerics unchanged (same arithmetic order).

"""
Extrapolate solution to one face of the reference element.

face ∈ (:west, :east, :south, :north) corresponding to ξ=-1, ξ=+1, η=-1, η=+1.
Returns matrix (Np, Neq) of face traces along the face SPs.
"""
function face_trace(u_e, ops::FROperators{T}, face::Symbol) where {T}
    Np = size(u_e, 1)
    Neq = size(u_e, 3)
    out = zeros(T, Np, Neq)
    face_trace!(out, u_e, ops, face)
    return out
end

"""In-place face trace into `out` of size (Np, Neq)."""
function face_trace!(out::AbstractMatrix{T}, u_e, ops::FROperators{T}, face::Symbol) where {T}
    Np = size(u_e, 1)
    Neq = size(u_e, 3)
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

@inline function _copy_comp!(dest::AbstractVector{T}, src, e::Int, q::Int, Neq::Int) where {T}
    @inbounds for c in 1:Neq
        dest[c] = src[q, c, e]
    end
    return dest
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
    ug::AbstractVector{T},
    fh::AbstractVector{T},
) where {T}
    ug_res = exterior_state(bc, u_int, nx_out, ny_out, x, y, t)
    @inbounds for c in 1:length(ug)
        ug[c] = ug_res[c]
    end
    if into_domain
        interface_flux_n!(fh, eq, ug, u_int, -nx_out, -ny_out, flux_kind)
    else
        interface_flux_n!(fh, eq, u_int, ug, nx_out, ny_out, flux_kind)
    end
    @inbounds for c in 1:length(fh)
        fhat_face[q, c, e] = fh[c] * sJ
    end
    return nothing
end

# Backward-compatible boundary_fhat! without work buffers
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
    Neq = length(u_int)
    return boundary_fhat!(
        fhat_face, q, e, eq, bc, u_int, nx_out, ny_out, sJ, x, y, t, flux_kind, into_domain,
        Vector{T}(undef, Neq), Vector{T}(undef, Neq),
    )
end

"""
    residual!(du, state::SolutionState2D, eq, method)

Metric-aware 2D strong-form FR residual (Cartesian or curved) with capturing hooks.
Supports Periodic / Transmissive / Reflecting / Dirichlet / GhostState BCs and
optional solid-element masks (forward-facing step etc.).

Phase 4: reuses work buffers within a residual evaluation; arithmetic order matches
the pre-optimization residual (no threading by default → bit-stable results).
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

    # Face solution traces (one allocation set per residual call)
    trW = zeros(T, Np, Neq, Nel)
    trE = zeros(T, Np, Neq, Nel)
    trS = zeros(T, Np, Neq, Nel)
    trN = zeros(T, Np, Neq, Nel)
    @inbounds for e in 1:Nel
        has_solid && is_solid(mesh, e) && continue
        u_e = @view u_work[:, :, e, :]
        face_trace!(@view(trW[:, :, e]), u_e, ops, :west)
        face_trace!(@view(trE[:, :, e]), u_e, ops, :east)
        face_trace!(@view(trS[:, :, e]), u_e, ops, :south)
        face_trace!(@view(trN[:, :, e]), u_e, ops, :north)
    end

    fhat_W = zeros(T, Np, Neq, Nel)
    fhat_E = zeros(T, Np, Neq, Nel)
    fhat_S = zeros(T, Np, Neq, Nel)
    fhat_N = zeros(T, Np, Neq, Nel)

    flux_kind = state.scheme.flux
    wall = ReflectingBC()

    # Scratch (reused for all faces / volume points)
    u_m = Vector{T}(undef, Neq)
    u_p = Vector{T}(undef, Neq)
    ug = Vector{T}(undef, Neq)
    fh = Vector{T}(undef, Neq)
    Fx = Vector{T}(undef, Neq)
    Gy = Vector{T}(undef, Neq)
    Ft = zeros(T, Np, Np, Neq)
    Gt = zeros(T, Np, Np, Neq)

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
                _copy_comp!(u_m, trW, eR, q, Neq)
                x, y = physical_xy(mesh, eR, -one(T), ops.ξ[q])
                boundary_fhat!(
                    fhat_W, q, eR, eq, wall, u_m,
                    met.nx_W[q, eR], met.ny_W[q, eR], met.sJ_W[q, eR],
                    x, y, t, flux_kind, true, ug, fh,
                )
            elseif sR && !sL
                _copy_comp!(u_m, trE, eL, q, Neq)
                x, y = physical_xy(mesh, eL, one(T), ops.ξ[q])
                boundary_fhat!(
                    fhat_E, q, eL, eq, wall, u_m, nnx, nny, sJ,
                    x, y, t, flux_kind, false, ug, fh,
                )
            else
                _copy_comp!(u_m, trE, eL, q, Neq)
                _copy_comp!(u_p, trW, eR, q, Neq)
                interface_flux_n!(fh, eq, u_m, u_p, nnx, nny, flux_kind)
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
                    _copy_comp!(u_m, trE, eR, q, Neq)
                    _copy_comp!(u_p, trW, eL, q, Neq)
                    nnx, nny, sJ = met.nx_E[q, eR], met.ny_E[q, eR], met.sJ_E[q, eR]
                    interface_flux_n!(fh, eq, u_m, u_p, nnx, nny, flux_kind)
                    for c in 1:Neq
                        val = fh[c] * sJ
                        fhat_E[q, c, eR] = val
                        fhat_W[q, c, eL] = val
                    end
                end
            else
                if !(has_solid && is_solid(mesh, eL))
                    _copy_comp!(u_m, trW, eL, q, Neq)
                    x, y = physical_xy(mesh, eL, -one(T), ops.ξ[q])
                    boundary_fhat!(
                        fhat_W, q, eL, eq, mesh.left_bc, u_m,
                        met.nx_W[q, eL], met.ny_W[q, eL], met.sJ_W[q, eL],
                        x, y, t, flux_kind, true, ug, fh,
                    )
                end
                if !(has_solid && is_solid(mesh, eR))
                    _copy_comp!(u_p, trE, eR, q, Neq)
                    xR, yR = physical_xy(mesh, eR, one(T), ops.ξ[q])
                    boundary_fhat!(
                        fhat_E, q, eR, eq, mesh.right_bc, u_p,
                        met.nx_E[q, eR], met.ny_E[q, eR], met.sJ_E[q, eR],
                        xR, yR, t, flux_kind, false, ug, fh,
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
                _copy_comp!(u_m, trS, eT, q, Neq)
                x, y = physical_xy(mesh, eT, ops.ξ[q], -one(T))
                boundary_fhat!(
                    fhat_S, q, eT, eq, wall, u_m,
                    met.nx_S[q, eT], met.ny_S[q, eT], met.sJ_S[q, eT],
                    x, y, t, flux_kind, true, ug, fh,
                )
            elseif sT && !sB
                _copy_comp!(u_m, trN, eB, q, Neq)
                x, y = physical_xy(mesh, eB, ops.ξ[q], one(T))
                boundary_fhat!(
                    fhat_N, q, eB, eq, wall, u_m, nnx, nny, sJ,
                    x, y, t, flux_kind, false, ug, fh,
                )
            else
                _copy_comp!(u_m, trN, eB, q, Neq)
                _copy_comp!(u_p, trS, eT, q, Neq)
                interface_flux_n!(fh, eq, u_m, u_p, nnx, nny, flux_kind)
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
                    _copy_comp!(u_m, trN, eT, q, Neq)
                    _copy_comp!(u_p, trS, eB, q, Neq)
                    nnx, nny, sJ = met.nx_N[q, eT], met.ny_N[q, eT], met.sJ_N[q, eT]
                    interface_flux_n!(fh, eq, u_m, u_p, nnx, nny, flux_kind)
                    for c in 1:Neq
                        val = fh[c] * sJ
                        fhat_N[q, c, eT] = val
                        fhat_S[q, c, eB] = val
                    end
                end
            else
                if !(has_solid && is_solid(mesh, eB))
                    _copy_comp!(u_m, trS, eB, q, Neq)
                    x, y = physical_xy(mesh, eB, ops.ξ[q], -one(T))
                    boundary_fhat!(
                        fhat_S, q, eB, eq, mesh.bottom_bc, u_m,
                        met.nx_S[q, eB], met.ny_S[q, eB], met.sJ_S[q, eB],
                        x, y, t, flux_kind, true, ug, fh,
                    )
                end
                if !(has_solid && is_solid(mesh, eT))
                    _copy_comp!(u_p, trN, eT, q, Neq)
                    xT, yT = physical_xy(mesh, eT, ops.ξ[q], one(T))
                    boundary_fhat!(
                        fhat_N, q, eT, eq, mesh.top_bc, u_p,
                        met.nx_N[q, eT], met.ny_N[q, eT], met.sJ_N[q, eT],
                        xT, yT, t, flux_kind, false, ug, fh,
                    )
                end
            end
        end
    end

    # Volume residual (reuse Ft, Gt; in-place physical fluxes)
    @inbounds for e in 1:Nel
        has_solid && is_solid(mesh, e) && continue
        for j in 1:Np, i in 1:Np
            Uij = @view u_work[i, j, e, :]
            physical_flux_x!(Fx, eq, Uij)
            physical_flux_y!(Gy, eq, Uij)
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
                fill!(@view(du[:, :, e, :]), zero(T))
            end
        end
    end
    return du
end

function residual!(du, state::SolutionState2D, eq)
    return residual!(du, state, eq, NullCapturing())
end
