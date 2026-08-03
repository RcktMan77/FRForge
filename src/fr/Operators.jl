# FROperators: reference-element discrete operators for 1D FR.

"""
    FROperators{T}

Reference-element operators for polynomial degree `p` on GL solution points.
"""
struct FROperators{T}
    p::Int
    ξ::Vector{T}       # GL SPs, length Np = p+1
    w::Vector{T}       # GL weights on [-1,1]
    D::Matrix{T}       # d/dξ at SPs
    ℓ_L::Vector{T}     # Lagrange basis at ξ=-1
    ℓ_R::Vector{T}     # Lagrange basis at ξ=+1
    gL_ξ::Vector{T}    # g'_DG,L(ξ_j)
    gR_ξ::Vector{T}    # g'_DG,R(ξ_j)
    gL::Vector{T}      # g_DG,L(ξ_j) (for tests / diagnostics)
    gR::Vector{T}      # g_DG,R(ξ_j)
end

"""
    build_operators(p; T=Float64) -> FROperators{T}

Build GL nodes/weights, differentiation matrix, Lagrange endpoint values,
and g_DG correction derivatives (Appendix C).
"""
function build_operators(p::Int; T::Type=Float64)
    p >= 0 || throw(ArgumentError("p must be >= 0"))
    Np = p + 1
    ξ, w = gauss_legendre_nodes_weights(Np; T=T)
    D = differentiation_matrix(ξ)
    ℓ_L = lagrange_at(ξ, -one(T))
    ℓ_R = lagrange_at(ξ, one(T))
    gL, gR, gL_ξ, gR_ξ = g_DG_values_and_derivs(p, ξ)
    return FROperators{T}(p, ξ, w, D, ℓ_L, ℓ_R, gL_ξ, gR_ξ, gL, gR)
end

"""Number of solution points per element."""
n_points(ops::FROperators) = length(ops.ξ)
