# Metric terms for 2D FR on (possibly curved) quadrilateral elements.
#
# Mapping x(ξ,η), y(ξ,η) on reference square [-1,1]².
# Jacobian J = x_ξ y_η - x_η y_ξ.
# Contravariant fluxes:
#   F̃ =  y_η F - x_η G
#   G̃ = -y_ξ F + x_ξ G
# so that ∇_x·(F,G) = (1/J)(∂_ξ F̃ + ∂_η G̃).

"""
    MeshMetrics2D{T}

Precomputed geometry/metrics at solution points and faces for a 2D mesh.
"""
struct MeshMetrics2D{T}
    x::Array{T,3}       # (Np, Np, Nel)
    y::Array{T,3}
    x_ξ::Array{T,3}
    x_η::Array{T,3}
    y_ξ::Array{T,3}
    y_η::Array{T,3}
    J::Array{T,3}       # det ∂(x,y)/∂(ξ,η)
    # Face: unit outward normal (nx,ny) and surface Jacobian sJ at face SPs
    # Vertical faces: index (q along η, e); Horizontal: (q along ξ, e)
    nx_W::Array{T,2}
    ny_W::Array{T,2}
    sJ_W::Array{T,2}
    nx_E::Array{T,2}
    ny_E::Array{T,2}
    sJ_E::Array{T,2}
    nx_S::Array{T,2}
    ny_S::Array{T,2}
    sJ_S::Array{T,2}
    nx_N::Array{T,2}
    ny_N::Array{T,2}
    sJ_N::Array{T,2}
    h_char::Vector{T}   # characteristic length per element (CFL / AV)
    is_curved::Bool
end

"""True if mesh carries high-order (non-affine) geometry."""
is_curved_mesh(mesh::Mesh2D) = mesh.geom_x !== nothing

"""
    physical_xy(mesh, e, ξ, η)

Physical coordinates. Uses high-order geometry when present, else affine.
"""
function physical_xy(mesh::Mesh2D{T}, e::Int, ξ::T, η::T) where {T}
    jx, jy = element_coords(mesh, e)
    xL, xR = mesh.x_vertices[jx], mesh.x_vertices[jx + 1]
    yB, yT = mesh.y_vertices[jy], mesh.y_vertices[jy + 1]
    X = (xL + xR) / T(2) + mesh.Jx[jx] * ξ
    Y = (yB + yT) / T(2) + mesh.Jy[jy] * η
    if mesh.wavy_amp !== nothing
        return wavy_physical(X, Y, mesh.wavy_amp)
    end
    if mesh.geom_x === nothing
        return X, Y
    end
    # Isoparametric Lagrange on geom nodes (same 1D node set as geom_ξ)
    ξg = mesh.geom_ξ
    Ng = length(ξg)
    ℓξ = lagrange_at(ξg, ξ)
    ℓη = lagrange_at(ξg, η)
    x = zero(T)
    y = zero(T)
    @inbounds for j in 1:Ng, i in 1:Ng
        w = ℓξ[i] * ℓη[j]
        x += w * mesh.geom_x[i, j, e]
        y += w * mesh.geom_y[i, j, e]
    end
    return x, y
end

"""
    build_mesh_metrics(mesh, ops) -> MeshMetrics2D

Build metrics at solution points. Affine Cartesian when `geom_x` is `nothing`;
isoparametric curved when geometry nodes are set.
"""
function build_mesh_metrics(mesh::Mesh2D{T}, ops::FROperators{T}) where {T}
    Np = n_points(ops)
    Nel = mesh.n_elements
    D = ops.D
    ξ = ops.ξ

    x = zeros(T, Np, Np, Nel)
    y = zeros(T, Np, Np, Nel)
    x_ξ = zeros(T, Np, Np, Nel)
    x_η = zeros(T, Np, Np, Nel)
    y_ξ = zeros(T, Np, Np, Nel)
    y_η = zeros(T, Np, Np, Nel)
    J = zeros(T, Np, Np, Nel)

    curved = mesh.geom_x !== nothing

    @inbounds for e in 1:Nel
        for j in 1:Np, i in 1:Np
            xx, yy = physical_xy(mesh, e, ξ[i], ξ[j])
            x[i, j, e] = xx
            y[i, j, e] = yy
        end
        # Differentiate mapping with FR differentiation matrix
        for j in 1:Np, i in 1:Np
            dxξ = zero(T)
            dyξ = zero(T)
            dxη = zero(T)
            dyη = zero(T)
            for k in 1:Np
                dxξ += D[i, k] * x[k, j, e]
                dyξ += D[i, k] * y[k, j, e]
                dxη += D[j, k] * x[i, k, e]
                dyη += D[j, k] * y[i, k, e]
            end
            x_ξ[i, j, e] = dxξ
            y_ξ[i, j, e] = dyξ
            x_η[i, j, e] = dxη
            y_η[i, j, e] = dyη
            J[i, j, e] = dxξ * dyη - dxη * dyξ
        end
    end

    # Face metrics via Lagrange extrapolation of volume metric fields
    nx_W = zeros(T, Np, Nel)
    ny_W = zeros(T, Np, Nel)
    sJ_W = zeros(T, Np, Nel)
    nx_E = zeros(T, Np, Nel)
    ny_E = zeros(T, Np, Nel)
    sJ_E = zeros(T, Np, Nel)
    nx_S = zeros(T, Np, Nel)
    ny_S = zeros(T, Np, Nel)
    sJ_S = zeros(T, Np, Nel)
    nx_N = zeros(T, Np, Nel)
    ny_N = zeros(T, Np, Nel)
    sJ_N = zeros(T, Np, Nel)
    ℓ_L, ℓ_R = ops.ℓ_L, ops.ℓ_R

    @inbounds for e in 1:Nel
        for j in 1:Np  # η-index on vertical faces
            # Extrapolate metric components to ξ=±1
            xη_W = zero(T)
            yη_W = zero(T)
            xη_E = zero(T)
            yη_E = zero(T)
            for i in 1:Np
                xη_W += ℓ_L[i] * x_η[i, j, e]
                yη_W += ℓ_L[i] * y_η[i, j, e]
                xη_E += ℓ_R[i] * x_η[i, j, e]
                yη_E += ℓ_R[i] * y_η[i, j, e]
            end
            # East (+ξ): n ∝ (y_η, -x_η), sJ = ||(y_η,-x_η)||
            sE = sqrt(yη_E * yη_E + xη_E * xη_E)
            sE = max(sE, eps(T))
            sJ_E[j, e] = sE
            nx_E[j, e] = yη_E / sE
            ny_E[j, e] = -xη_E / sE
            # West (-ξ): outward opposite
            sW = sqrt(yη_W * yη_W + xη_W * xη_W)
            sW = max(sW, eps(T))
            sJ_W[j, e] = sW
            nx_W[j, e] = -yη_W / sW
            ny_W[j, e] = xη_W / sW
        end
        for i in 1:Np  # ξ-index on horizontal faces
            xξ_S = zero(T)
            yξ_S = zero(T)
            xξ_N = zero(T)
            yξ_N = zero(T)
            for j in 1:Np
                xξ_S += ℓ_L[j] * x_ξ[i, j, e]
                yξ_S += ℓ_L[j] * y_ξ[i, j, e]
                xξ_N += ℓ_R[j] * x_ξ[i, j, e]
                yξ_N += ℓ_R[j] * y_ξ[i, j, e]
            end
            # North (+η): n ∝ (-y_ξ, x_ξ)
            sN = sqrt(yξ_N * yξ_N + xξ_N * xξ_N)
            sN = max(sN, eps(T))
            sJ_N[i, e] = sN
            nx_N[i, e] = -yξ_N / sN
            ny_N[i, e] = xξ_N / sN
            # South (-η)
            sS = sqrt(yξ_S * yξ_S + xξ_S * xξ_S)
            sS = max(sS, eps(T))
            sJ_S[i, e] = sS
            nx_S[i, e] = yξ_S / sS
            ny_S[i, e] = -xξ_S / sS
        end
    end

    h_char = zeros(T, Nel)
    @inbounds for e in 1:Nel
        # Characteristic length ~ min edge length ≈ 2 * min(sJ) / (something)
        # Use sqrt of mean |J| * 4 / max(1,p) as area-based scale, and min face sJ
        area = zero(T)
        for j in 1:Np, i in 1:Np
            area += ops.w[i] * ops.w[j] * abs(J[i, j, e])
        end
        h_area = sqrt(max(area, eps(T)))
        h_face = typemax(T)
        for q in 1:Np
            h_face = min(h_face, sJ_W[q, e], sJ_E[q, e], sJ_S[q, e], sJ_N[q, e])
        end
        # sJ ~ physical edge length / 2 for reference length 2; h ~ sJ
        h_char[e] = min(h_area, T(2) * h_face)
        if !curved
            jx, jy = element_coords(mesh, e)
            h_char[e] = min(mesh.Δx[jx], mesh.Δy[jy])
        end
    end

    return MeshMetrics2D{T}(
        x, y, x_ξ, x_η, y_ξ, y_η, J,
        nx_W, ny_W, sJ_W, nx_E, ny_E, sJ_E,
        nx_S, ny_S, sJ_S, nx_N, ny_N, sJ_N,
        h_char, curved,
    )
end

"""
    apply_geometry_warp!(mesh, ops, warp!)

Set isoparametric geometry nodes at solution-point locations by warping
the affine image: `(x,y) = warp(x_aff, y_aff)`.
"""
function apply_geometry_warp!(
    mesh::Mesh2D{T},
    ops::FROperators{T},
    warp!::Function,
) where {T}
    Np = n_points(ops)
    Nel = mesh.n_elements
    gx = zeros(T, Np, Np, Nel)
    gy = zeros(T, Np, Np, Nel)
    @inbounds for e in 1:Nel
        for j in 1:Np, i in 1:Np
            # Affine position (ignore geom)
            jx, jy = element_coords(mesh, e)
            xL, xR = mesh.x_vertices[jx], mesh.x_vertices[jx + 1]
            yB, yT = mesh.y_vertices[jy], mesh.y_vertices[jy + 1]
            xa = (xL + xR) / T(2) + mesh.Jx[jx] * ops.ξ[i]
            ya = (yB + yT) / T(2) + mesh.Jy[jy] * ops.ξ[j]
            xw, yw = warp!(xa, ya)
            gx[i, j, e] = xw
            gy[i, j, e] = yw
        end
    end
    mesh.geom_x = gx
    mesh.geom_y = gy
    mesh.geom_ξ = copy(ops.ξ)
    return mesh
end

"""
    make_wavy_mesh2d(nx, ny, ops; amp=0.05, kwargs...) -> Mesh2D

Unit-square structured mesh with sinusoidal warp (classic free-stream / order test).
Geometry nodes are isoparametric; use `build_mesh_metrics_analytic_wavy` for
free-stream-preserving analytic metrics (recommended for residual).
"""
function make_wavy_mesh2d(
    nx::Int,
    ny::Int,
    ops::FROperators{T};
    amp::Real=0.05,
    x_left=0.0,
    x_right=1.0,
    y_bottom=0.0,
    y_top=1.0,
    left_bc::AbstractBC=PeriodicBC(),
    right_bc::AbstractBC=PeriodicBC(),
    bottom_bc::AbstractBC=PeriodicBC(),
    top_bc::AbstractBC=PeriodicBC(),
) where {T}
    mesh = Mesh2D(
        x_left,
        x_right,
        y_bottom,
        y_top,
        nx,
        ny;
        left_bc=left_bc,
        right_bc=right_bc,
        bottom_bc=bottom_bc,
        top_bc=top_bc,
        T=T,
    )
    a = T(amp)
    apply_geometry_warp!(mesh, ops, (x, y) -> wavy_physical(x, y, a))
    mesh.wavy_amp = a
    return mesh
end

"""
    wavy_physical(X, Y, amp) -> (x, y)

Global wavy map from computational (X,Y) to physical (x,y).
"""
function wavy_physical(X::T, Y::T, amp::T) where {T}
    x = X + amp * sin(T(2π) * Y) * sin(T(π) * X)
    y = Y + amp * sin(T(2π) * X) * sin(T(π) * Y)
    return x, y
end

"""
    wavy_partials(X, Y, amp) -> (x_X, x_Y, y_X, y_Y)

Analytic Jacobian of the wavy map (computational → physical).
"""
function wavy_partials(X::T, Y::T, amp::T) where {T}
    # x = X + a sin(2π Y) sin(π X)
    x_X = one(T) + amp * sin(T(2π) * Y) * T(π) * cos(T(π) * X)
    x_Y = amp * T(2π) * cos(T(2π) * Y) * sin(T(π) * X)
    # y = Y + a sin(2π X) sin(π Y)
    y_X = amp * T(2π) * cos(T(2π) * X) * sin(T(π) * Y)
    y_Y = one(T) + amp * sin(T(2π) * X) * T(π) * cos(T(π) * Y)
    return x_X, x_Y, y_X, y_Y
end

"""
    build_mesh_metrics_analytic_wavy(mesh, ops; amp) -> MeshMetrics2D

Analytic metrics via chain rule through affine computational coords.
Face metrics are continuous across structured interfaces → free-stream
preservation holds (GCL). Prefer this over pure discrete D-metrics for wavy tests.
"""
function build_mesh_metrics_analytic_wavy(
    mesh::Mesh2D{T},
    ops::FROperators{T};
    amp::Real=0.05,
) where {T}
    Np = n_points(ops)
    Nel = mesh.n_elements
    a = T(amp)
    ξ = ops.ξ
    D = ops.D

    x = zeros(T, Np, Np, Nel)
    y = zeros(T, Np, Np, Nel)
    x_ξ = zeros(T, Np, Np, Nel)
    x_η = zeros(T, Np, Np, Nel)
    y_ξ = zeros(T, Np, Np, Nel)
    y_η = zeros(T, Np, Np, Nel)
    J = zeros(T, Np, Np, Nel)

    @inbounds for e in 1:Nel
        jx, jy = element_coords(mesh, e)
        Jx, Jy = mesh.Jx[jx], mesh.Jy[jy]
        xL, xR = mesh.x_vertices[jx], mesh.x_vertices[jx + 1]
        yB, yT = mesh.y_vertices[jy], mesh.y_vertices[jy + 1]
        for j in 1:Np, i in 1:Np
            X = (xL + xR) / T(2) + Jx * ξ[i]
            Y = (yB + yT) / T(2) + Jy * ξ[j]
            xx, yy = wavy_physical(X, Y, a)
            x[i, j, e] = xx
            y[i, j, e] = yy
            x_X, x_Y, y_X, y_Y = wavy_partials(X, Y, a)
            # chain rule: ∂/∂ξ = ∂/∂X * Jx, ∂/∂η = ∂/∂Y * Jy
            x_ξ[i, j, e] = x_X * Jx
            x_η[i, j, e] = x_Y * Jy
            y_ξ[i, j, e] = y_X * Jx
            y_η[i, j, e] = y_Y * Jy
            J[i, j, e] = x_ξ[i, j, e] * y_η[i, j, e] - x_η[i, j, e] * y_ξ[i, j, e]
        end
    end

    # Face metrics: evaluate analytic map on faces (continuous across elements)
    nx_W = zeros(T, Np, Nel)
    ny_W = zeros(T, Np, Nel)
    sJ_W = zeros(T, Np, Nel)
    nx_E = zeros(T, Np, Nel)
    ny_E = zeros(T, Np, Nel)
    sJ_E = zeros(T, Np, Nel)
    nx_S = zeros(T, Np, Nel)
    ny_S = zeros(T, Np, Nel)
    sJ_S = zeros(T, Np, Nel)
    nx_N = zeros(T, Np, Nel)
    ny_N = zeros(T, Np, Nel)
    sJ_N = zeros(T, Np, Nel)

    @inbounds for e in 1:Nel
        jx, jy = element_coords(mesh, e)
        Jx, Jy = mesh.Jx[jx], mesh.Jy[jy]
        xL, xR = mesh.x_vertices[jx], mesh.x_vertices[jx + 1]
        yB, yT = mesh.y_vertices[jy], mesh.y_vertices[jy + 1]
        for j in 1:Np
            Y = (yB + yT) / T(2) + Jy * ξ[j]
            # West ξ=-1
            Xw = (xL + xR) / T(2) + Jx * (-one(T))
            x_X, x_Y, y_X, y_Y = wavy_partials(Xw, Y, a)
            xη = x_Y * Jy
            yη = y_Y * Jy
            sW = max(sqrt(yη * yη + xη * xη), eps(T))
            sJ_W[j, e] = sW
            nx_W[j, e] = -yη / sW
            ny_W[j, e] = xη / sW
            # East ξ=+1
            Xe = (xL + xR) / T(2) + Jx * (one(T))
            x_X, x_Y, y_X, y_Y = wavy_partials(Xe, Y, a)
            xη = x_Y * Jy
            yη = y_Y * Jy
            sE = max(sqrt(yη * yη + xη * xη), eps(T))
            sJ_E[j, e] = sE
            nx_E[j, e] = yη / sE
            ny_E[j, e] = -xη / sE
        end
        for i in 1:Np
            X = (xL + xR) / T(2) + Jx * ξ[i]
            # South η=-1
            Ys = (yB + yT) / T(2) + Jy * (-one(T))
            x_X, x_Y, y_X, y_Y = wavy_partials(X, Ys, a)
            xξ = x_X * Jx
            yξ = y_X * Jx
            sS = max(sqrt(yξ * yξ + xξ * xξ), eps(T))
            sJ_S[i, e] = sS
            nx_S[i, e] = yξ / sS
            ny_S[i, e] = -xξ / sS
            # North η=+1
            Yn = (yB + yT) / T(2) + Jy * (one(T))
            x_X, x_Y, y_X, y_Y = wavy_partials(X, Yn, a)
            xξ = x_X * Jx
            yξ = y_X * Jx
            sN = max(sqrt(yξ * yξ + xξ * xξ), eps(T))
            sJ_N[i, e] = sN
            nx_N[i, e] = -yξ / sN
            ny_N[i, e] = xξ / sN
        end
    end

    h_char = zeros(T, Nel)
    @inbounds for e in 1:Nel
        jx, jy = element_coords(mesh, e)
        h_char[e] = min(mesh.Δx[jx], mesh.Δy[jy])
    end

    return MeshMetrics2D{T}(
        x, y, x_ξ, x_η, y_ξ, y_η, J,
        nx_W, ny_W, sJ_W, nx_E, ny_E, sJ_E,
        nx_S, ny_S, sJ_S, nx_N, ny_N, sJ_N,
        h_char, true,
    )
end
