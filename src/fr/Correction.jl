# g_DG left/right correction functions (Legendre/Radau) + Legendre Vandermonde.

"""
    legendre_P(k, ξ) -> P_k(ξ)

Standard Legendre polynomial P_k on [-1,1] via three-term recurrence.
"""
function legendre_P(k::Int, ξ::T) where {T<:Number}
    k < 0 && throw(ArgumentError("k must be >= 0"))
    if k == 0
        return one(T)
    elseif k == 1
        return ξ
    end
    Pkm2 = one(T)
    Pkm1 = ξ
    Pk = Pkm1
    for m in 2:k
        Pk = ((T(2m - 1) * ξ * Pkm1) - (T(m - 1) * Pkm2)) / T(m)
        Pkm2 = Pkm1
        Pkm1 = Pk
    end
    return Pk
end

"""
    legendre_P_prime(k, ξ) -> P_k'(ξ)

Analytic derivative. Endpoints: P_k'(1)=k(k+1)/2, P_k'(-1)=(-1)^{k+1} k(k+1)/2.
Interior: (1-ξ²) P_k' = k P_{k-1} - k ξ P_k.
"""
function legendre_P_prime(k::Int, ξ::T) where {T<:Number}
    k < 0 && throw(ArgumentError("k must be >= 0"))
    k == 0 && return zero(T)
    # Endpoint formulas (stable)
    if ξ >= one(T) - 1000 * eps(float(one(T)))
        return T(k * (k + 1)) / T(2)
    elseif ξ <= -one(T) + 1000 * eps(float(one(T)))
        # (-1)^{k+1} * k(k+1)/2
        return (isodd(k) ? one(T) : -one(T)) * T(k * (k + 1)) / T(2)
    end
    Pk = legendre_P(k, ξ)
    Pkm1 = legendre_P(k - 1, ξ)
    return (T(k) * Pkm1 - T(k) * ξ * Pk) / (one(T) - ξ * ξ)
end

"""
    g_DG_values_and_derivs(p, ξ) -> (gL, gR, gL_ξ, gR_ξ)

Construct g_L, g_R and their derivatives at solution points ξ:

  r_R = P_{p+1} - P_p   (vanishes at +1)
  r_L = P_{p+1} + P_p   (vanishes at -1)
  g_L = r_R / r_R(-1)
  g_R = r_L / r_L(+1)
"""
function g_DG_values_and_derivs(p::Int, ξ::AbstractVector{T}) where {T}
    p >= 0 || throw(ArgumentError("p must be >= 0"))
    n = length(ξ)

    rR_m1 = legendre_P(p + 1, -one(T)) - legendre_P(p, -one(T))
    rL_p1 = legendre_P(p + 1, one(T)) + legendre_P(p, one(T))
    abs(rR_m1) < 100 * eps(T) && error("r_R(-1) vanished for p=$p")
    abs(rL_p1) < 100 * eps(T) && error("r_L(+1) vanished for p=$p")

    gL = zeros(T, n)
    gR = zeros(T, n)
    gL_ξ = zeros(T, n)
    gR_ξ = zeros(T, n)

    for (i, x) in enumerate(ξ)
        Pp = legendre_P(p, x)
        Pp1 = legendre_P(p + 1, x)
        rR = Pp1 - Pp
        rL = Pp1 + Pp
        gL[i] = rR / rR_m1
        gR[i] = rL / rL_p1

        dPp = legendre_P_prime(p, x)
        dPp1 = legendre_P_prime(p + 1, x)
        gL_ξ[i] = (dPp1 - dPp) / rR_m1
        gR_ξ[i] = (dPp1 + dPp) / rL_p1
    end

    return gL, gR, gL_ξ, gR_ξ
end

"""
    legendre_vandermonde(ξ) -> V

V[j,k] = P_{k-1}(ξ_j) for k=1..Np (modal basis for Persson sensor).
"""
function legendre_vandermonde(ξ::AbstractVector{T}) where {T}
    Np = length(ξ)
    V = zeros(T, Np, Np)
    @inbounds for j in 1:Np, k in 1:Np
        V[j, k] = legendre_P(k - 1, ξ[j])
    end
    return V
end

"""
    g_DG_endpoints(p; T=Float64) -> NamedTuple

Evaluate g_L, g_R at ξ=±1 for unit tests.
"""
function g_DG_endpoints(p::Int; T::Type = Float64)
    ξ = T[-one(T), one(T)]
    gL, gR, _, _ = g_DG_values_and_derivs(p, ξ)
    return (
        gL_m1 = gL[1],
        gL_p1 = gL[2],
        gR_m1 = gR[1],
        gR_p1 = gR[2],
    )
end
