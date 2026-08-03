# Cartesian / curved 2D mesh of quadrilateral elements.

"""
    Mesh2D{T}

Structured `nx × ny` quad mesh. Affine Cartesian by default; optional
isoparametric high-order geometry via `geom_x`, `geom_y`, `geom_ξ`.

Element index `e = (jy - 1) * nx + jx` with `jx ∈ 1:nx`, `jy ∈ 1:ny`.
"""
mutable struct Mesh2D{T}
    nx::Int
    ny::Int
    n_elements::Int
    x_vertices::Vector{T}   # length nx+1 (computational / affine skeleton)
    y_vertices::Vector{T}   # length ny+1
    Δx::Vector{T}           # length nx
    Δy::Vector{T}           # length ny
    Jx::Vector{T}           # Δx/2 per column
    Jy::Vector{T}           # Δy/2 per row
    left_bc::AbstractBC
    right_bc::AbstractBC
    bottom_bc::AbstractBC
    top_bc::AbstractBC
    # Isoparametric geometry (nothing = pure affine from vertices)
    geom_x::Union{Nothing,Array{T,3}}  # (Ng, Ng, Nel)
    geom_y::Union{Nothing,Array{T,3}}
    geom_ξ::Union{Nothing,Vector{T}}   # 1D nodes for geom Lagrange
    # If set, physical_xy uses analytic wavy map with this amplitude
    wavy_amp::Union{Nothing,T}
end

function Mesh2D(
    x_left,
    x_right,
    y_bottom,
    y_top,
    nx::Int,
    ny::Int;
    left_bc::AbstractBC=PeriodicBC(),
    right_bc::AbstractBC=PeriodicBC(),
    bottom_bc::AbstractBC=PeriodicBC(),
    top_bc::AbstractBC=PeriodicBC(),
    T::Type=Float64,
)
    nx >= 1 && ny >= 1 || throw(ArgumentError("nx, ny must be >= 1"))
    x_right > x_left || throw(ArgumentError("x_right > x_left required"))
    y_top > y_bottom || throw(ArgumentError("y_top > y_bottom required"))
    if left_bc isa PeriodicBC || right_bc isa PeriodicBC
        (left_bc isa PeriodicBC && right_bc isa PeriodicBC) ||
            throw(ArgumentError("x-periodic requires both left and right PeriodicBC"))
    end
    if bottom_bc isa PeriodicBC || top_bc isa PeriodicBC
        (bottom_bc isa PeriodicBC && top_bc isa PeriodicBC) ||
            throw(ArgumentError("y-periodic requires both bottom and top PeriodicBC"))
    end
    xv = collect(range(T(x_left), T(x_right); length=nx + 1))
    yv = collect(range(T(y_bottom), T(y_top); length=ny + 1))
    Δx = diff(xv)
    Δy = diff(yv)
    return Mesh2D{T}(
        nx,
        ny,
        nx * ny,
        xv,
        yv,
        Δx,
        Δy,
        Δx ./ T(2),
        Δy ./ T(2),
        left_bc,
        right_bc,
        bottom_bc,
        top_bc,
        nothing,
        nothing,
        nothing,
        nothing,
    )
end

"""Element index from (jx, jy) both 1-based."""
@inline element_index(mesh::Mesh2D, jx::Int, jy::Int) = (jy - 1) * mesh.nx + jx

"""(jx, jy) from linear element index e (1-based)."""
function element_coords(mesh::Mesh2D, e::Int)
    jy, r = divrem(e - 1, mesh.nx)
    return r + 1, jy + 1
end

"""All SP physical coordinates as (X, Y) arrays size (Np, Np, Nel)."""
function physical_coords_2d(mesh::Mesh2D{T}, ops::FROperators{T}) where {T}
    Np = n_points(ops)
    Nel = mesh.n_elements
    X = zeros(T, Np, Np, Nel)
    Y = zeros(T, Np, Np, Nel)
    for e in 1:Nel
        for jy in 1:Np, jx in 1:Np
            x, y = physical_xy(mesh, e, ops.ξ[jx], ops.ξ[jy])
            X[jx, jy, e] = x
            Y[jx, jy, e] = y
        end
    end
    return X, Y
end
