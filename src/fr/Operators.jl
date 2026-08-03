# FROperators: reference-element discrete operators for 1D FR.

"""
    gauss_lobatto_legendre_nodes_weights(n; T=Float64) -> (ξ, w)

`n` Gauss–Lobatto–Legendre nodes and weights on `[-1, 1]`.
Nodes are endpoints ±1 plus interior roots of ``P'_p`` where ``p = n-1``.
Weights: ``w_i = 2 / (p(p+1) [P_p(ξ_i)]²)`` (then normalized to sum 2).
"""
function gauss_lobatto_legendre_nodes_weights(n::Int; T::Type=Float64)
    n >= 1 || throw(ArgumentError("n must be >= 1, got $n"))
    if n == 1
        return T[0], T[2]
    end
    p = n - 1
    ξ = Vector{T}(undef, n)
    ξ[1] = -one(T)
    ξ[n] = one(T)
    # Interior: Newton on P'_p = 0, Chebyshev–Lobatto initial guess
    for j in 2:(n - 1)
        x = cos(T(π) * T(j - 1) / T(p))
        for _ in 1:80
            Pp = legendre_P(p, x)
            dPp = legendre_P_prime(p, x)
            # (1-x²) P'' = 2x P' - p(p+1) P
            one_m_x2 = one(T) - x * x
            if abs(one_m_x2) < 1000 * eps(T)
                break
            end
            d2Pp = (T(2) * x * dPp - T(p * (p + 1)) * Pp) / one_m_x2
            if abs(d2Pp) < eps(T)
                break
            end
            xnew = x - dPp / d2Pp
            xnew = min(max(xnew, -one(T) + 10 * eps(T)), one(T) - 10 * eps(T))
            if abs(xnew - x) < 10 * eps(T)
                x = xnew
                break
            end
            x = xnew
        end
        ξ[j] = x
    end
    sort!(ξ)
    ξ[1] = -one(T)
    ξ[n] = one(T)
    w = Vector{T}(undef, n)
    for i in 1:n
        Pp = legendre_P(p, ξ[i])
        den = T(p) * T(p + 1) * Pp * Pp
        w[i] = T(2) / max(den, eps(T))
    end
    w .*= T(2) / sum(w)
    return ξ, w
end

"""
    FROperators{T}

Reference-element operators for polynomial degree `p` on GL or GLL solution points.
"""
struct FROperators{T}
    p::Int
    points::Symbol     # :gl | :gll
    ξ::Vector{T}       # SPs, length Np = p+1
    w::Vector{T}       # quadrature weights on [-1,1]
    D::Matrix{T}       # d/dξ at SPs
    ℓ_L::Vector{T}     # Lagrange basis at ξ=-1
    ℓ_R::Vector{T}     # Lagrange basis at ξ=+1
    gL_ξ::Vector{T}    # g'_DG,L(ξ_j)
    gR_ξ::Vector{T}    # g'_DG,R(ξ_j)
    gL::Vector{T}      # g_DG,L(ξ_j) (for tests / diagnostics)
    gR::Vector{T}      # g_DG,R(ξ_j)
end

"""
    build_operators(p; points=:gl, T=Float64) -> FROperators{T}

Build nodes/weights, differentiation matrix, Lagrange endpoint values,
and g_DG correction derivatives (Appendix C).

`points`: `:gl` (Gauss–Legendre, default) or `:gll` (Gauss–Lobatto–Legendre).
"""
function build_operators(p::Int; points::Symbol=:gl, T::Type=Float64)
    p >= 0 || throw(ArgumentError("p must be >= 0"))
    points in (:gl, :gll) || throw(ArgumentError("points must be :gl or :gll, got $points"))
    Np = p + 1
    if points === :gl
        ξ, w = gauss_legendre_nodes_weights(Np; T=T)
    else
        ξ, w = gauss_lobatto_legendre_nodes_weights(Np; T=T)
    end
    D = differentiation_matrix(ξ)
    ℓ_L = lagrange_at(ξ, -one(T))
    ℓ_R = lagrange_at(ξ, one(T))
    gL, gR, gL_ξ, gR_ξ = g_DG_values_and_derivs(p, ξ)
    return FROperators{T}(p, points, ξ, w, D, ℓ_L, ℓ_R, gL_ξ, gR_ξ, gL, gR)
end

build_operators(p::Int, scheme::SchemeConfig; T::Type=Float64) =
    build_operators(p; points=scheme.points, T=T)

"""Number of solution points per element."""
n_points(ops::FROperators) = length(ops.ξ)
