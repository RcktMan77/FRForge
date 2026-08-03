# Capturing / discontinuity-treatment hook pipeline (null defaults for M1).

abstract type AbstractCapturingMethod end
abstract type AbstractShockSensor end
abstract type AbstractDissipationOperator end

struct NullSensor <: AbstractShockSensor end
struct NullDissipation <: AbstractDissipationOperator end

"""Default method: identity hooks. Used from M1 onward."""
struct NullCapturing <: AbstractCapturingMethod
    sensor::NullSensor
    dissip::NullDissipation
end
NullCapturing() = NullCapturing(NullSensor(), NullDissipation())

"""Traces at element interfaces: uL[e,c], uR[e,c]."""
struct InterfaceTraces{T}
    uL::Matrix{T}  # (Nel, Neq)
    uR::Matrix{T}
end

function allocate_traces(Nel::Int, Neq::Int, ::Type{T}=Float64) where {T}
    return InterfaceTraces{T}(zeros(T, Nel, Neq), zeros(T, Nel, Neq))
end

"""Optional copy/limit of conserved state before residual."""
function preprocess_state!(u_work, method::AbstractCapturingMethod, state, eq)
    copyto!(u_work, state.u)
    return u_work
end

"""Build left/right traces per element via Lagrange ℓ(±1)."""
function extrapolate_interface!(traces::InterfaceTraces, method::AbstractCapturingMethod, u_work, state, eq)
    ops = state.ops
    Np, Nel, Neq = size(u_work)
    @inbounds for e in 1:Nel
        for c in 1:Neq
            sL = zero(eltype(u_work))
            sR = zero(eltype(u_work))
            for j in 1:Np
                sL += ops.ℓ_L[j] * u_work[j, e, c]
                sR += ops.ℓ_R[j] * u_work[j, e, c]
            end
            traces.uL[e, c] = sL
            traces.uR[e, c] = sR
        end
    end
    return traces
end

"""Optional numerical flux override; return `nothing` to use equation default."""
function numerical_flux_method(method::AbstractCapturingMethod, eq, uL, uR)
    return nothing
end

function sense!(σ, method::AbstractCapturingMethod, u_work, state, eq)
    fill!(σ, zero(eltype(σ)))
    return σ
end

function apply_dissipation!(du, method::AbstractCapturingMethod, σ, u_work, state, eq)
    return du
end

function post_step!(state, method::AbstractCapturingMethod, eq)
    return nothing
end
