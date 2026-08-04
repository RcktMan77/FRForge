# Exact Riemann solution for the Sod shock tube (ideal gas).

"""
Classic Sod left/right states.
  Left:  ρ=1, u=0, p=1
  Right: ρ=0.125, u=0, p=0.1
  γ=1.4, interface at x0 (default 0.5), typically t=0.2 on [0,1]
"""
struct SodProblem{T}
    ρL::T
    uL::T
    pL::T
    ρR::T
    uR::T
    pR::T
    γ::T
    x0::T
end

function SodProblem(; γ = 1.4, x0 = 0.5, T::Type = Float64)
    return SodProblem{T}(T(1), T(0), T(1), T(0.125), T(0), T(0.1), T(γ), T(x0))
end

_cs(γ, ρ, p) = sqrt(γ * p / ρ)

function sod_p_star(prob::SodProblem{T}; tol = T(1e-12), maxiter = 60) where {T}
    γ = prob.γ
    ρL, uL, pL = prob.ρL, prob.uL, prob.pL
    ρR, uR, pR = prob.ρR, prob.uR, prob.pR
    cL = _cs(γ, ρL, pL)
    cR = _cs(γ, ρR, pR)
    p = T(0.5) * (pL + pR)
    p = max(p, T(1e-8))

    function f_df(p)
        if p > pL
            AL = T(2) / ((γ + 1) * ρL)
            BL = (γ - 1) / (γ + 1) * pL
            ql = sqrt(AL / (p + BL))
            fL = (p - pL) * ql
            dfL = ql * (one(T) - T(0.5) * (p - pL) / (p + BL))
        else
            fL = 2cL / (γ - 1) * ((p / pL)^((γ - 1) / (2γ)) - one(T))
            dfL = (1 / (ρL * cL)) * (p / pL)^(-(γ + 1) / (2γ))
        end
        if p > pR
            AR = T(2) / ((γ + 1) * ρR)
            BR = (γ - 1) / (γ + 1) * pR
            qr = sqrt(AR / (p + BR))
            fR = (p - pR) * qr
            dfR = qr * (one(T) - T(0.5) * (p - pR) / (p + BR))
        else
            fR = 2cR / (γ - 1) * ((p / pR)^((γ - 1) / (2γ)) - one(T))
            dfR = (1 / (ρR * cR)) * (p / pR)^(-(γ + 1) / (2γ))
        end
        return fL + fR + (uR - uL), dfL + dfR
    end

    for _ in 1:maxiter
        f, df = f_df(p)
        abs(f) < tol && break
        p = max(p - f / df, tol)
    end
    return p
end

function sod_u_star(prob::SodProblem{T}, p_star) where {T}
    γ = prob.γ
    ρL, uL, pL = prob.ρL, prob.uL, prob.pL
    ρR, uR, pR = prob.ρR, prob.uR, prob.pR
    cL = _cs(γ, ρL, pL)
    cR = _cs(γ, ρR, pR)
    if p_star > pL
        AL = T(2) / ((γ + 1) * ρL)
        BL = (γ - 1) / (γ + 1) * pL
        fL = (p_star - pL) * sqrt(AL / (p_star + BL))
    else
        fL = 2cL / (γ - 1) * ((p_star / pL)^((γ - 1) / (2γ)) - one(T))
    end
    if p_star > pR
        AR = T(2) / ((γ + 1) * ρR)
        BR = (γ - 1) / (γ + 1) * pR
        fR = (p_star - pR) * sqrt(AR / (p_star + BR))
    else
        fR = 2cR / (γ - 1) * ((p_star / pR)^((γ - 1) / (2γ)) - one(T))
    end
    return T(0.5) * (uL + uR) + T(0.5) * (fR - fL)
end

function _rho_star(γ, ρ, p, p_star)
    if p_star > p
        return ρ * (p_star / p + (γ - 1) / (γ + 1)) / ((γ - 1) / (γ + 1) * p_star / p + 1)
    else
        return ρ * (p_star / p)^(1 / γ)
    end
end

"""
    sod_exact(prob, x, t) -> (ρ, u, p)

Exact primitives at (x, t) for the Sod problem.
"""
function sod_exact(prob::SodProblem{T}, x::Real, t::Real) where {T}
    tT = T(t)
    xT = T(x)
    if tT <= zero(T)
        return xT < prob.x0 ? (prob.ρL, prob.uL, prob.pL) : (prob.ρR, prob.uR, prob.pR)
    end
    γ = prob.γ
    ρL, uL, pL = prob.ρL, prob.uL, prob.pL
    ρR, uR, pR = prob.ρR, prob.uR, prob.pR
    cL = _cs(γ, ρL, pL)
    cR = _cs(γ, ρR, pR)
    p_s = sod_p_star(prob)
    u_s = sod_u_star(prob, p_s)
    ξ = (xT - prob.x0) / tT
    ρ_sL = _rho_star(γ, ρL, pL, p_s)
    ρ_sR = _rho_star(γ, ρR, pR, p_s)

    # Left side of contact
    if ξ < u_s
        if p_s > pL  # left shock
            sL = uL - cL * sqrt((γ + 1) / (2γ) * p_s / pL + (γ - 1) / (2γ))
            return ξ < sL ? (ρL, uL, pL) : (ρ_sL, u_s, p_s)
        else  # left rarefaction
            sHL = uL - cL
            c_sL = cL * (p_s / pL)^((γ - 1) / (2γ))
            sTL = u_s - c_sL
            if ξ < sHL
                return ρL, uL, pL
            elseif ξ < sTL
                u = (2 / (γ + 1)) * (cL + ((γ - 1) / 2) * uL + ξ)
                c = (2 / (γ + 1)) * (cL + ((γ - 1) / 2) * (uL - ξ))
                ρ = ρL * (c / cL)^(2 / (γ - 1))
                p = pL * (ρ / ρL)^γ
                return ρ, u, p
            else
                return ρ_sL, u_s, p_s
            end
        end
    end

    # Right of contact
    if p_s > pR  # right shock
        sR = uR + cR * sqrt((γ + 1) / (2γ) * p_s / pR + (γ - 1) / (2γ))
        return ξ > sR ? (ρR, uR, pR) : (ρ_sR, u_s, p_s)
    else
        sHR = uR + cR
        c_sR = cR * (p_s / pR)^((γ - 1) / (2γ))
        sTR = u_s + c_sR
        if ξ > sHR
            return ρR, uR, pR
        elseif ξ > sTR
            u = (2 / (γ + 1)) * (-cR + ((γ - 1) / 2) * uR + ξ)
            c = (2 / (γ + 1)) * (cR - ((γ - 1) / 2) * (uR - ξ))
            ρ = ρR * (c / cR)^(2 / (γ - 1))
            p = pR * (ρ / ρR)^γ
            return ρ, u, p
        else
            return ρ_sR, u_s, p_s
        end
    end
end

function sod_exact_conserved(eq::Euler1D, prob::SodProblem, x, t)
    ρ, u, p = sod_exact(prob, x, t)
    return primitives_to_conserved(eq, ρ, u, p)
end
