# Persson–Peraire modal sensor + element-local artificial viscosity (invent baseline).

"""
    PerssonSensor{T}

Modal smoothness sensor (Persson–Peraire). Operates on density for Euler and
on the scalar field for Burgers/advection.
"""
struct PerssonSensor{T} <: AbstractShockSensor
    κ::T
    s0_factor::T   # s0 = s0_factor * log10(p); default -4
    ε_floor::T
end

function PerssonSensor(;
    κ::Real = 4.0,
    s0_factor::Real = -4.0,
    ε_floor::Real = 1e-16,
    T::Type = Float64,
)
    return PerssonSensor{T}(T(κ), T(s0_factor), T(ε_floor))
end

"""
    ElementArtificialViscosity{T}

Element-scalar artificial viscosity:
  ε_e = c_av * σ_e * (h_e / p) * λ_max,e

`av_form`:
- `"conservative_br0"` — FR lifting of viscous flux g=ε∂_x u with BR0 averages + IP (default)
- `"element_local_DD"` — r += (ε/J²) D(Du)  (legacy; may not conserve with variable ε)
"""
struct ElementArtificialViscosity{T} <: AbstractDissipationOperator
    c_av::T
    av_form::String
end

function ElementArtificialViscosity(;
    c_av::Real = 0.1,
    av_form::AbstractString = "conservative_br0",
    T::Type = Float64,
)
    return ElementArtificialViscosity{T}(T(c_av), String(av_form))
end

"""
    PerssonAVMethod{T}

Baseline capturing method: Persson modal sensor + element AV.
Residual never names this type; only hooks are used.
"""
struct PerssonAVMethod{T} <: AbstractCapturingMethod
    sensor::PerssonSensor{T}
    dissip::ElementArtificialViscosity{T}
end

function PerssonAVMethod(;
    κ::Real = 4.0,
    s0_factor::Real = -4.0,
    c_av::Real = 0.1,
    av_form::AbstractString = "conservative_br0",
    T::Type = Float64,
)
    return PerssonAVMethod{T}(
        PerssonSensor(; κ = κ, s0_factor = s0_factor, T = T),
        ElementArtificialViscosity(; c_av = c_av, av_form = av_form, T = T),
    )
end

# c_av=0.1 is the explicit-RK-friendly scored default (0.5 can stiffen BR0).
default_persson_params() = (
    κ = 4.0,
    s0_factor = -4.0,
    c_av = 0.1,
    av_form = "conservative_br0",
)

function method_params(m::PerssonAVMethod)
    return Dict{String,Any}(
        "κ" => m.sensor.κ,
        "s0_factor" => m.sensor.s0_factor,
        "c_av" => m.dissip.c_av,
        "av_form" => m.dissip.av_form,
        "method" => "persson_av",
    )
end

method_params(::NullCapturing) = Dict{String,Any}("method" => "null")
method_params(::AbstractCapturingMethod) = Dict{String,Any}()

# --- Method-level hook forwarding (composition, not field probing in residual) ---

function sense!(σ, method::PerssonAVMethod, u_work, state, eq)
    return sense!(σ, method.sensor, u_work, state, eq)
end

function apply_dissipation!(du, method::PerssonAVMethod, σ, u_work, state, eq)
    return apply_dissipation!(du, method.dissip, σ, u_work, state, eq)
end

"""Sensor scalar field: density for Euler, first component otherwise."""
function sensor_field(::Euler1D, u_work::AbstractArray{T,3}, e::Int) where {T}
    return @view u_work[:, e, 1]
end

function sensor_field(::AbstractEquation, u_work::AbstractArray{T,3}, e::Int) where {T}
    return @view u_work[:, e, 1]
end

"""
    sense!(σ, sensor::PerssonSensor, u_work, state, eq)

Fill element-wise σ ∈ (0,1) from modal energy indicator.
"""
function sense!(σ::AbstractVector{T}, sensor::PerssonSensor{T}, u_work, state, eq) where {T}
    ops = state.ops
    p = ops.p
    Np = length(ops.ξ)
    Nel = size(u_work, 2)
    length(σ) == Nel || throw(DimensionMismatch("σ length $(length(σ)) != Nel $Nel"))

    V = ops.V_legendre
    # Modal coeffs û = V \\ u  (same algorithm; V cached on operators)
    s0 = sensor.s0_factor * log10(T(max(p, 1)))

    @inbounds for e in 1:Nel
        u_e = collect(sensor_field(eq, u_work, e))
        û = V \ u_e
        e_high = û[Np]^2
        e_tot = sum(abs2, û) + sensor.ε_floor
        s_e = log10(e_high / e_tot + sensor.ε_floor)
        # Smooth indicator
        σ[e] = one(T) / (one(T) + exp(-sensor.κ * (s_e - s0)))
    end
    return σ
end

"""Per-element max wave speed for viscosity scaling."""
function element_max_wavespeed(eq, u_work::AbstractArray{T,3}, e::Int) where {T}
    Np = size(u_work, 1)
    m = zero(T)
    if eq isa Euler1D
        @inbounds for j in 1:Np
            U = @view u_work[j, e, :]
            m = max(m, abs(velocity(eq, U)) + sound_speed(eq, U))
        end
    else
        @inbounds for j in 1:Np
            m = max(m, abs(u_work[j, e, 1]))
        end
        if eq isa LinearAdvection1D
            m = max(m, abs(eq.a))
        end
    end
    return max(m, eps(T))
end

"""Compute element viscosities ε_e from σ and local wave speeds."""
function element_viscosities(
    dissip::ElementArtificialViscosity{T},
    σ::AbstractVector{T},
    u_work,
    state,
    eq,
) where {T}
    mesh, ops = state.mesh, state.ops
    p = max(ops.p, 1)
    Nel = size(u_work, 2)
    ε = zeros(T, Nel)
    @inbounds for e in 1:Nel
        σe = σ[e]
        σe <= zero(T) && continue
        h = mesh.Δx[e]
        λ = element_max_wavespeed(eq, u_work, e)
        ε[e] = dissip.c_av * σe * (h / T(p)) * λ
    end
    return ε
end

"""
    apply_dissipation!(du, dissip::ElementArtificialViscosity, σ, u_work, state, eq)

Add AV residual. Default `conservative_br0` uses FR lifting of viscous flux
g = ε ∂_x u with BR0 averages and interior penalty so discrete mass telescopes
on periodic meshes.
"""
function apply_dissipation!(
    du::AbstractArray{T,3},
    dissip::ElementArtificialViscosity{T},
    σ::AbstractVector{T},
    u_work,
    state,
    eq,
) where {T}
    if dissip.av_form == "conservative_br0"
        return apply_dissipation_br0!(du, dissip, σ, u_work, state, eq)
    elseif dissip.av_form == "element_local_DD"
        return apply_dissipation_local_DD!(du, dissip, σ, u_work, state, eq)
    else
        error("Unsupported av_form=$(dissip.av_form)")
    end
end

function apply_dissipation_local_DD!(
    du::AbstractArray{T,3},
    dissip::ElementArtificialViscosity{T},
    σ::AbstractVector{T},
    u_work,
    state,
    eq,
) where {T}
    mesh, ops = state.mesh, state.ops
    p = max(ops.p, 1)
    Np, Nel, Neq = size(u_work)
    D = ops.D
    ε = element_viscosities(dissip, σ, u_work, state, eq)
    v = zeros(T, Np)
    @inbounds for e in 1:Nel
        εe = ε[e]
        εe <= zero(T) && continue
        J = mesh.J[e]
        scale = εe / (J * J)
        for c in 1:Neq
            for j in 1:Np
                s = zero(T)
                for k in 1:Np
                    s += D[j, k] * u_work[k, e, c]
                end
                v[j] = s
            end
            for j in 1:Np
                s = zero(T)
                for k in 1:Np
                    s += D[j, k] * v[k]
                end
                du[j, e, c] += scale * s
            end
        end
    end
    return du
end

"""
Conservative BR0 viscous residual for u_t += ∂_x(ε ∂_x u):

  g = ε (1/J) D u   (discontinuous viscous flux at SPs)
  interface: ĝ = ½(g⁻+g⁺) − τ (u⁺−u⁻),  τ ∼ ε̄ (p+1)² / h̄
  residual:  du += (1/J)[ D g + (ĝ_L−g_L) g'_L + (ĝ_R−g_R) g'_R ]
"""
function apply_dissipation_br0!(
    du::AbstractArray{T,3},
    dissip::ElementArtificialViscosity{T},
    σ::AbstractVector{T},
    u_work,
    state,
    eq,
) where {T}
    mesh, ops = state.mesh, state.ops
    Np, Nel, Neq = size(u_work)
    D, ℓ_L, ℓ_R = ops.D, ops.ℓ_L, ops.ℓ_R
    gL_ξ, gR_ξ = ops.gL_ξ, ops.gR_ξ
    p = ops.p
    rws = ensure_residual_workspace!(state)
    b = ensure_br0_workspace_1d!(rws, Np, Nel, Neq)
    ε = b.ε
    # fill viscosities into pooled buffer (same values as element_viscosities)
    ε_new = element_viscosities(dissip, σ, u_work, state, eq)
    copyto!(ε, ε_new)

    # g[j,e,c] viscous flux at SPs; also interface traces of u and g
    g, uL, uR, gL, gR = b.g, b.uL, b.uR, b.gL, b.gR
    ghat_L, ghat_R = b.ghat_L, b.ghat_R

    @inbounds for e in 1:Nel
        J = mesh.J[e]
        εe = ε[e]
        for c in 1:Neq
            # q = (1/J) D u; g = ε q
            for j in 1:Np
                dq = zero(T)
                for k in 1:Np
                    dq += D[j, k] * u_work[k, e, c]
                end
                g[j, e, c] = εe * dq / J
            end
            s_uL = zero(T)
            s_uR = zero(T)
            s_gL = zero(T)
            s_gR = zero(T)
            for j in 1:Np
                s_uL += ℓ_L[j] * u_work[j, e, c]
                s_uR += ℓ_R[j] * u_work[j, e, c]
                s_gL += ℓ_L[j] * g[j, e, c]
                s_gR += ℓ_R[j] * g[j, e, c]
            end
            uL[e, c] = s_uL
            uR[e, c] = s_uR
            gL[e, c] = s_gL
            gR[e, c] = s_gR
        end
    end

    function br0_flux(e_left::Int, e_right::Int, c::Int)
        # states from left element right face / right element left face
        um = uR[e_left, c]
        up = uL[e_right, c]
        gm = gR[e_left, c]
        gp = gL[e_right, c]
        εm = ε[e_left]
        εp = ε[e_right]
        hm = mesh.Δx[e_left]
        hp = mesh.Δx[e_right]
        εbar = T(0.5) * (εm + εp)
        hbar = T(0.5) * (hm + hp)
        # Moderate interior penalty (full (p+1)²/h can over-stiffen explicit RK)
        τ = T(2) * εbar * T(p + 1) / max(hbar, eps(T))
        return T(0.5) * (gm + gp) - τ * (up - um)
    end

    @inbounds for e in 1:(Nel - 1)
        for c in 1:Neq
            gh = br0_flux(e, e + 1, c)
            ghat_R[e, c] = gh
            ghat_L[e + 1, c] = gh
        end
    end

    # Domain ends: periodic or Neumann-like (zero viscous flux jump / freestream)
    if mesh.left_bc isa PeriodicBC
        @inbounds for c in 1:Neq
            gh = br0_flux(Nel, 1, c)
            ghat_R[Nel, c] = gh
            ghat_L[1, c] = gh
        end
    else
        # Transmissive/Dirichlet: zero viscous numerical flux (no penalty) using interior
        @inbounds for c in 1:Neq
            ghat_L[1, c] = gL[1, c]
            ghat_R[Nel, c] = gR[Nel, c]
        end
    end

    # FR strong form for +∂_x g  (added to hyperbolic residual which already has -∂_x F)
    @inbounds for e in 1:Nel
        J = mesh.J[e]
        for c in 1:Neq
            for j in 1:Np
                vol = zero(T)
                for k in 1:Np
                    vol += D[j, k] * g[k, e, c]
                end
                corr =
                    (ghat_L[e, c] - gL[e, c]) * gL_ξ[j] +
                    (ghat_R[e, c] - gR[e, c]) * gR_ξ[j]
                du[j, e, c] += (vol + corr) / J
            end
        end
    end
    return du
end

"""
Mass residual of viscous operator on a constant field (should be ~0).
Uses one residual evaluation with σ≡1, constant u.
"""
function viscous_mass_residual_scale(ops::FROperators{T}, mesh::Mesh1D{T}; p = ops.p) where {T}
    # Constant state → gradients zero → BR0 residual zero
    state = allocate_state(mesh, ops, Val(1))
    fill!(state.u, one(T))
    eq = Burgers1D()
    dissip = ElementArtificialViscosity(; c_av = 1.0, av_form = "conservative_br0", T = T)
    σ = ones(T, mesh.n_elements)
    du = zeros(T, size(state.u)...)
    apply_dissipation!(du, dissip, σ, state.u, state, eq)
    mass = zero(T)
    Np, Nel = size(state.u, 1), size(state.u, 2)
    for e in 1:Nel, j in 1:Np
        mass += ops.w[j] * mesh.J[e] * du[j, e, 1]
    end
    return abs(mass)
end
