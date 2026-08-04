# 2D Persson modal sensor + element artificial viscosity (Cartesian / curved quads).
# Used by invent baseline and docs; residual never names PerssonAVMethod.

"""Sensor scalar field on element `e`: density for Euler2D, first component otherwise."""
function sensor_field_2d(::Euler2D, u_work::AbstractArray{T,4}, e::Int) where {T}
    return @view u_work[:, :, e, 1]
end

function sensor_field_2d(::AbstractEquation, u_work::AbstractArray{T,4}, e::Int) where {T}
    return @view u_work[:, :, e, 1]
end

"""
    sense!(σ, sensor::PerssonSensor, u_work::Array{T,4}, state, eq)

2D Persson–Peraire modal sensor on tensor-product elements.
Modal coefficients ``\\hat{U} = V^{-1} U V^{-T}``; high-mode energy uses the
highest Legendre mode in either ξ or η.
"""
function sense!(
    σ::AbstractVector{T},
    sensor::PerssonSensor{T},
    u_work::AbstractArray{T,4},
    state::SolutionState2D,
    eq,
) where {T}
    ops = state.ops
    p = ops.p
    Np = length(ops.ξ)
    Nel = size(u_work, 3)
    length(σ) == Nel || throw(DimensionMismatch("σ length $(length(σ)) != Nel $Nel"))

    V = ops.V_legendre
    s0 = sensor.s0_factor * log10(T(max(p, 1)))

    # Element-local sensor (read-only cached V); parallel when FRFORGE residual threads > 1
    foreach_element(Nel) do e
        U = Matrix{T}(sensor_field_2d(eq, u_work, e))  # (Np, Np)
        # û = V^{-1} U V^{-T}
        Û = V \ U
        Û = Û / V'
        e_tot = sum(abs2, Û) + sensor.ε_floor
        e_high = zero(T)
        @inbounds for j in 1:Np
            e_high += Û[Np, j]^2
        end
        @inbounds for i in 1:(Np - 1)
            e_high += Û[i, Np]^2
        end
        s_e = log10(e_high / e_tot + sensor.ε_floor)
        σ[e] = one(T) / (one(T) + exp(-sensor.κ * (s_e - s0)))
    end
    return σ
end

"""Per-element max wave speed for 2D viscosity scaling."""
function element_max_wavespeed_2d(eq, u_work::AbstractArray{T,4}, e::Int) where {T}
    Np = size(u_work, 1)
    m = zero(T)
    if eq isa Euler2D
        @inbounds for j in 1:Np, i in 1:Np
            U = @view u_work[i, j, e, :]
            ρ = max(U[1], eps(T))
            u = U[2] / ρ
            v = U[3] / ρ
            c = sound_speed(eq, U)
            m = max(m, abs(u) + c, abs(v) + c)
        end
    elseif eq isa LinearAdvection2D
        m = hypot(eq.ax, eq.ay)
    else
        @inbounds for j in 1:Np, i in 1:Np
            m = max(m, abs(u_work[i, j, e, 1]))
        end
    end
    return max(m, eps(T))
end

"""
Element viscosities on Mesh2D (Cartesian or curved):
  ε_e = c_av * σ_e * (h_e / p) * λ_max,e
with h_e from metric characteristic length (`metrics.h_char`).
"""
function element_viscosities_2d(
    dissip::ElementArtificialViscosity{T},
    σ::AbstractVector{T},
    u_work::AbstractArray{T,4},
    state::SolutionState2D,
    eq,
) where {T}
    ops = state.ops
    p = max(ops.p, 1)
    Nel = size(u_work, 3)
    ε = zeros(T, Nel)
    hchar = state.metrics.h_char
    @inbounds for e in 1:Nel
        σe = σ[e]
        σe <= zero(T) && continue
        h = hchar[e]
        λ = element_max_wavespeed_2d(eq, u_work, e)
        ε[e] = dissip.c_av * σe * (h / T(p)) * λ
    end
    return ε
end

"""
    apply_dissipation!(du::Array{T,4}, dissip, σ, u_work, state, eq)

2D AV. `conservative_br0` uses tensor-product BR0 in ξ and η (Cartesian).
`element_local_DD` uses element-local Laplacian in reference space.
"""
function apply_dissipation!(
    du::AbstractArray{T,4},
    dissip::ElementArtificialViscosity{T},
    σ::AbstractVector{T},
    u_work::AbstractArray{T,4},
    state::SolutionState2D,
    eq,
) where {T}
    if dissip.av_form == "conservative_br0"
        return apply_dissipation_br0_2d!(du, dissip, σ, u_work, state, eq)
    elseif dissip.av_form == "element_local_DD"
        return apply_dissipation_local_DD_2d!(du, dissip, σ, u_work, state, eq)
    else
        error("Unsupported av_form=$(dissip.av_form) for 2D")
    end
end

"""Element-local Laplacian AV: du += ε (∂²u/∂x² + ∂²u/∂y²) via reference D²."""
function apply_dissipation_local_DD_2d!(
    du::AbstractArray{T,4},
    dissip::ElementArtificialViscosity{T},
    σ::AbstractVector{T},
    u_work::AbstractArray{T,4},
    state::SolutionState2D,
    eq,
) where {T}
    mesh, ops = state.mesh, state.ops
    Np = size(u_work, 1)
    Nel = size(u_work, 3)
    Neq = size(u_work, 4)
    D = ops.D
    ε = element_viscosities_2d(dissip, σ, u_work, state, eq)
    v_pool = make_vector_scratch_pool(T, Np)
    foreach_element(Nel) do e
        εe = ε[e]
        εe <= zero(T) && return
        v = vector_scratch_for(v_pool)
        jx, jy = element_coords(mesh, e)
        Jx, Jy = mesh.Jx[jx], mesh.Jy[jy]
        sx = εe / (Jx * Jx)
        sy = εe / (Jy * Jy)
        @inbounds for c in 1:Neq
            for j in 1:Np
                for i in 1:Np
                    s = zero(T)
                    for k in 1:Np
                        s += D[i, k] * u_work[k, j, e, c]
                    end
                    v[i] = s
                end
                for i in 1:Np
                    s = zero(T)
                    for k in 1:Np
                        s += D[i, k] * v[k]
                    end
                    du[i, j, e, c] += sx * s
                end
            end
            for i in 1:Np
                for j in 1:Np
                    s = zero(T)
                    for k in 1:Np
                        s += D[j, k] * u_work[i, k, e, c]
                    end
                    v[j] = s
                end
                for j in 1:Np
                    s = zero(T)
                    for k in 1:Np
                        s += D[j, k] * v[k]
                    end
                    du[i, j, e, c] += sy * s
                end
            end
        end
    end
    return du
end

"""
2D BR0-style AV on Cartesian quads (tensor product):

  g^x = ε (1/Jx) ∂_ξ u ,   g^y = ε (1/Jy) ∂_η u
  interface averages + IP on vertical (x-normal) and horizontal (y-normal) faces
  FR lift for +∇·g added to hyperbolic residual.
"""
function apply_dissipation_br0_2d!(
    du::AbstractArray{T,4},
    dissip::ElementArtificialViscosity{T},
    σ::AbstractVector{T},
    u_work::AbstractArray{T,4},
    state::SolutionState2D,
    eq,
) where {T}
    mesh, ops = state.mesh, state.ops
    Np = size(u_work, 1)
    Nel = size(u_work, 3)
    Neq = size(u_work, 4)
    nx, ny = mesh.nx, mesh.ny
    D, ℓ_L, ℓ_R = ops.D, ops.ℓ_L, ops.ℓ_R
    gL_ξ, gR_ξ = ops.gL_ξ, ops.gR_ξ
    p = ops.p
    ε = element_viscosities_2d(dissip, σ, u_work, state, eq)

    rws = ensure_residual_workspace!(state)
    b = ensure_br0_workspace!(rws, Np, Nel, Neq)
    gx, gy = b.gx, b.gy
    uW, uE, uS, uN = b.uW, b.uE, b.uS, b.uN
    gxW, gxE, gyS, gyN = b.gxW, b.gxE, b.gyS, b.gyN
    ghat_W, ghat_E, ghat_S, ghat_N = b.ghat_W, b.ghat_E, b.ghat_S, b.ghat_N

    # Element-local gx/gy + face traces (parallel-safe writes)
    foreach_element(Nel) do e
        jx, jy = element_coords(mesh, e)
        Jx, Jy = mesh.Jx[jx], mesh.Jy[jy]
        εe = ε[e]
        @inbounds for c in 1:Neq
            for j in 1:Np, i in 1:Np
                dξ = zero(T)
                dη = zero(T)
                for k in 1:Np
                    dξ += D[i, k] * u_work[k, j, e, c]
                    dη += D[j, k] * u_work[i, k, e, c]
                end
                gx[i, j, e, c] = εe * dξ / Jx
                gy[i, j, e, c] = εe * dη / Jy
            end
            for j in 1:Np
                suW = zero(T)
                suE = zero(T)
                sgW = zero(T)
                sgE = zero(T)
                for i in 1:Np
                    suW += ℓ_L[i] * u_work[i, j, e, c]
                    suE += ℓ_R[i] * u_work[i, j, e, c]
                    sgW += ℓ_L[i] * gx[i, j, e, c]
                    sgE += ℓ_R[i] * gx[i, j, e, c]
                end
                uW[j, e, c] = suW
                uE[j, e, c] = suE
                gxW[j, e, c] = sgW
                gxE[j, e, c] = sgE
            end
            for i in 1:Np
                suS = zero(T)
                suN = zero(T)
                sgS = zero(T)
                sgN = zero(T)
                for j in 1:Np
                    suS += ℓ_L[j] * u_work[i, j, e, c]
                    suN += ℓ_R[j] * u_work[i, j, e, c]
                    sgS += ℓ_L[j] * gy[i, j, e, c]
                    sgN += ℓ_R[j] * gy[i, j, e, c]
                end
                uS[i, e, c] = suS
                uN[i, e, c] = suN
                gyS[i, e, c] = sgS
                gyN[i, e, c] = sgN
            end
        end
    end

    function br0_x(eL::Int, eR::Int, q::Int, c::Int)
        um = uE[q, eL, c]
        up = uW[q, eR, c]
        gm = gxE[q, eL, c]
        gp = gxW[q, eR, c]
        jxL, _ = element_coords(mesh, eL)
        jxR, _ = element_coords(mesh, eR)
        εbar = T(0.5) * (ε[eL] + ε[eR])
        hbar = T(0.5) * (mesh.Δx[jxL] + mesh.Δx[jxR])
        τ = T(2) * εbar * T(p + 1) / max(hbar, eps(T))
        return T(0.5) * (gm + gp) - τ * (up - um)
    end

    function br0_y(eB::Int, eT::Int, q::Int, c::Int)
        um = uN[q, eB, c]
        up = uS[q, eT, c]
        gm = gyN[q, eB, c]
        gp = gyS[q, eT, c]
        _, jyB = element_coords(mesh, eB)
        _, jyT = element_coords(mesh, eT)
        εbar = T(0.5) * (ε[eB] + ε[eT])
        hbar = T(0.5) * (mesh.Δy[jyB] + mesh.Δy[jyT])
        τ = T(2) * εbar * T(p + 1) / max(hbar, eps(T))
        return T(0.5) * (gm + gp) - τ * (up - um)
    end

    # Interior vertical faces
    @inbounds for jy in 1:ny, jx in 1:(nx - 1)
        eL = element_index(mesh, jx, jy)
        eR = element_index(mesh, jx + 1, jy)
        for q in 1:Np, c in 1:Neq
            gh = br0_x(eL, eR, q, c)
            ghat_E[q, eL, c] = gh
            ghat_W[q, eR, c] = gh
        end
    end
    # Domain left/right
    @inbounds for jy in 1:ny
        eL = element_index(mesh, 1, jy)
        eR = element_index(mesh, nx, jy)
        for q in 1:Np, c in 1:Neq
            if mesh.left_bc isa PeriodicBC
                gh = br0_x(eR, eL, q, c)
                ghat_E[q, eR, c] = gh
                ghat_W[q, eL, c] = gh
            else
                ghat_W[q, eL, c] = gxW[q, eL, c]
                ghat_E[q, eR, c] = gxE[q, eR, c]
            end
        end
    end
    # Interior horizontal faces
    @inbounds for jy in 1:(ny - 1), jx in 1:nx
        eB = element_index(mesh, jx, jy)
        eT = element_index(mesh, jx, jy + 1)
        for q in 1:Np, c in 1:Neq
            gh = br0_y(eB, eT, q, c)
            ghat_N[q, eB, c] = gh
            ghat_S[q, eT, c] = gh
        end
    end
    # Domain bottom/top
    @inbounds for jx in 1:nx
        eB = element_index(mesh, jx, 1)
        eT = element_index(mesh, jx, ny)
        for q in 1:Np, c in 1:Neq
            if mesh.bottom_bc isa PeriodicBC
                gh = br0_y(eT, eB, q, c)
                ghat_N[q, eT, c] = gh
                ghat_S[q, eB, c] = gh
            else
                ghat_S[q, eB, c] = gyS[q, eB, c]
                ghat_N[q, eT, c] = gyN[q, eT, c]
            end
        end
    end

    # FR lift: du += ∇·g  (x and y contributions)
    @inbounds for e in 1:Nel
        jx, jy = element_coords(mesh, e)
        Jx, Jy = mesh.Jx[jx], mesh.Jy[jy]
        for c in 1:Neq
            # x-direction FR for each fixed j
            for j in 1:Np
                gL_face = zero(T)
                gR_face = zero(T)
                for i in 1:Np
                    gL_face += ℓ_L[i] * gx[i, j, e, c]
                    gR_face += ℓ_R[i] * gx[i, j, e, c]
                end
                for i in 1:Np
                    vol = zero(T)
                    for k in 1:Np
                        vol += D[i, k] * gx[k, j, e, c]
                    end
                    corr =
                        (ghat_W[j, e, c] - gL_face) * gL_ξ[i] +
                        (ghat_E[j, e, c] - gR_face) * gR_ξ[i]
                    du[i, j, e, c] += (vol + corr) / Jx
                end
            end
            # y-direction FR for each fixed i
            for i in 1:Np
                gS_face = zero(T)
                gN_face = zero(T)
                for j in 1:Np
                    gS_face += ℓ_L[j] * gy[i, j, e, c]
                    gN_face += ℓ_R[j] * gy[i, j, e, c]
                end
                for j in 1:Np
                    vol = zero(T)
                    for k in 1:Np
                        vol += D[j, k] * gy[i, k, e, c]
                    end
                    corr =
                        (ghat_S[i, e, c] - gS_face) * gL_ξ[j] +
                        (ghat_N[i, e, c] - gN_face) * gR_ξ[j]
                    du[i, j, e, c] += (vol + corr) / Jy
                end
            end
        end
    end
    return du
end

"""Mass residual of 2D viscous operator on a constant field (should be ~0)."""
function viscous_mass_residual_scale_2d(
    ops::FROperators{T},
    mesh::Mesh2D{T};
    p=ops.p,
) where {T}
    state = allocate_state(mesh, ops, Val(1))
    fill!(state.u, one(T))
    eq = LinearAdvection2D(one(T), zero(T))
    dissip = ElementArtificialViscosity(; c_av=1.0, av_form="conservative_br0", T=T)
    σ = ones(T, mesh.n_elements)
    du = zeros(T, size(state.u)...)
    apply_dissipation!(du, dissip, σ, state.u, state, eq)
    mass = zero(T)
    Np = n_points(ops)
    for e in 1:mesh.n_elements
        jx, jy = element_coords(mesh, e)
        Jx, Jy = mesh.Jx[jx], mesh.Jy[jy]
        for j in 1:Np, i in 1:Np
            mass += ops.w[i] * ops.w[j] * Jx * Jy * du[i, j, e, 1]
        end
    end
    return abs(mass)
end
