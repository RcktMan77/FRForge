# Uniform / nonuniform 1D mesh with BC tags.

"""
    Mesh1D{T}

1D mesh: vertices, element sizes, Jacobians J = Δx/2, left/right BCs.
"""
struct Mesh1D{T}
    n_elements::Int
    x_vertices::Vector{T}
    Δx::Vector{T}
    J::Vector{T}
    left_bc::AbstractBC
    right_bc::AbstractBC
end

"""
    Mesh1D(x_left, x_right, n_elements; left_bc, right_bc)

Uniform mesh on [x_left, x_right] with `n_elements` elements.
If either BC is periodic, both must be `PeriodicBC`.
"""
function Mesh1D(
    x_left,
    x_right,
    n_elements::Int;
    left_bc::AbstractBC = PeriodicBC(),
    right_bc::AbstractBC = PeriodicBC(),
    T::Type = Float64,
)
    n_elements >= 1 || throw(ArgumentError("n_elements must be >= 1"))
    x_right > x_left || throw(ArgumentError("x_right must be > x_left"))
    if is_periodic(left_bc) || is_periodic(right_bc)
        (left_bc isa PeriodicBC && right_bc isa PeriodicBC) ||
            throw(ArgumentError("Periodic BC must be set on both ends"))
    end

    xv = collect(range(T(x_left), T(x_right); length = n_elements + 1))
    Δx = diff(xv)
    J = Δx ./ T(2)
    return Mesh1D{T}(n_elements, xv, Δx, J, left_bc, right_bc)
end

"""Physical solution-point coordinates for element `e` (1-based)."""
function physical_coords(mesh::Mesh1D{T}, ops::FROperators{T}, e::Int) where {T}
    xL = mesh.x_vertices[e]
    xR = mesh.x_vertices[e + 1]
    # ξ ∈ [-1,1] → x = (xL+xR)/2 + (xR-xL)/2 * ξ
    mid = (xL + xR) / T(2)
    half = mesh.J[e]
    return mid .+ half .* ops.ξ
end

"""All physical SP coordinates as matrix (Np, Nel)."""
function physical_coords(mesh::Mesh1D{T}, ops::FROperators{T}) where {T}
    Np = n_points(ops)
    Nel = mesh.n_elements
    X = zeros(T, Np, Nel)
    for e in 1:Nel
        X[:, e] .= physical_coords(mesh, ops, e)
    end
    return X
end
