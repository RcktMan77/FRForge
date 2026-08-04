# FR reference points & Lagrange tools (1D).
# Layout (fr/): Points → Correction → Operators → Threading → ResidualWorkspace → Residual(1D/2D).

"""
    gauss_legendre_nodes_weights(n; T=Float64) -> (ξ, w)

Return `n` Gauss–Legendre nodes and weights on the reference interval `[-1, 1]`.
Uses the Golub–Welsch eigenvalue algorithm (symmetric Jacobi matrix).
"""
function gauss_legendre_nodes_weights(n::Int; T::Type=Float64)
    n >= 1 || throw(ArgumentError("n must be >= 1, got $n"))
    if n == 1
        return T[0], T[2]
    end

    # Jacobi matrix for Legendre: β_k = k / sqrt(4k² - 1)
    β = [T(k) / sqrt(T(4 * k * k - 1)) for k in 1:(n - 1)]
    J = zeros(T, n, n)
    for k in 1:(n - 1)
        J[k, k + 1] = β[k]
        J[k + 1, k] = β[k]
    end
    F = eigen(Symmetric(J))
    ξ = F.values
    # Sort ascending
    perm = sortperm(ξ)
    ξ = ξ[perm]
    V = F.vectors[:, perm]
    # Weights: w_i = 2 v_{1,i}^2  (first component of normalized eigenvector)
    w = T[2 * V[1, i]^2 for i in 1:n]
    # Normalize weights to sum to 2 (interval length)
    w .*= T(2) / sum(w)
    return ξ, w
end

"""
    lagrange_basis_matrix(ξ, x_eval) -> L

L[i,j] = ℓ_j(x_eval[i]) where ℓ_j are Lagrange basis polynomials on nodes ξ.
"""
function lagrange_basis_matrix(ξ::AbstractVector{T}, x_eval::AbstractVector{T}) where {T}
    n = length(ξ)
    m = length(x_eval)
    L = zeros(T, m, n)
    for i in 1:m
        x = x_eval[i]
        for j in 1:n
            ℓ = one(T)
            for k in 1:n
                if k != j
                    ℓ *= (x - ξ[k]) / (ξ[j] - ξ[k])
                end
            end
            L[i, j] = ℓ
        end
    end
    return L
end

"""
    lagrange_at(ξ, x) -> ℓ

Lagrange basis values ℓ_j(x) on nodes ξ, length n.
"""
function lagrange_at(ξ::AbstractVector{T}, x::T) where {T}
    n = length(ξ)
    ℓ = ones(T, n)
    for j in 1:n
        for k in 1:n
            if k != j
                ℓ[j] *= (x - ξ[k]) / (ξ[j] - ξ[k])
            end
        end
    end
    return ℓ
end

"""
    differentiation_matrix(ξ) -> D

D[i,j] = ℓ_j'(ξ_i) for Lagrange basis on nodes ξ.
"""
function differentiation_matrix(ξ::AbstractVector{T}) where {T}
    n = length(ξ)
    D = zeros(T, n, n)
    for i in 1:n
        for j in 1:n
            if i != j
                # barycentric-style product form for ℓ_j'(ξ_i)
                num = one(T)
                den = ξ[j] - ξ[i]
                for k in 1:n
                    if k != i && k != j
                        num *= (ξ[i] - ξ[k]) / (ξ[j] - ξ[k])
                    end
                end
                D[i, j] = num / den
            end
        end
        # Diagonal from partition of unity: sum_j ℓ_j' = 0
        D[i, i] = -sum(D[i, k] for k in 1:n if k != i)
    end
    return D
end
