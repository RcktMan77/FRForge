# Capturing / discontinuity-treatment hook pipeline.
#
# Residual stages (fixed order in residual!):
#   1. preprocess_state!        — positivity limiters, scaling, etc.
#   2. extrapolate_interface!   — hybrid FR–WENO / limiters on traces
#   3. numerical_flux_method    — optional flux override (else equation default)
#   4. FR volume + DG correction (fixed core; not replaceable via hooks)
#   5. sense! + apply_dissipation!  — AV / residual filters
#   6. post_step! (after each full time step, not RK substage) — solution limiters
#
# residual! dispatches only on AbstractCapturingMethod / these hooks —
# never on concrete method type names (e.g. never PerssonAVMethod by name).

abstract type AbstractCapturingMethod end
abstract type AbstractShockSensor end
abstract type AbstractDissipationOperator end

struct NullSensor <: AbstractShockSensor end
struct NullDissipation <: AbstractDissipationOperator end

"""Default method: identity hooks (no sensor or artificial viscosity)."""
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

# ---------------------------------------------------------------------------
# Hook API — defaults are true no-ops / standard FR behavior
# ---------------------------------------------------------------------------

"""Optional copy/limit of conserved state before residual."""
function preprocess_state!(u_work, method::AbstractCapturingMethod, state, eq)
    copyto!(u_work, state.u)
    return u_work
end

"""Build left/right traces per element via Lagrange ℓ(±1). Override for hybrid reconstruction."""
function extrapolate_interface!(
    traces::InterfaceTraces,
    method::AbstractCapturingMethod,
    u_work,
    state,
    eq,
)
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

"""Method-level sensor hook (default: zeros). Monolithic methods override this."""
function sense!(σ, method::AbstractCapturingMethod, u_work, state, eq)
    fill!(σ, zero(eltype(σ)))
    return σ
end

"""Method-level dissipation hook (default: no-op)."""
function apply_dissipation!(du, method::AbstractCapturingMethod, σ, u_work, state, eq)
    return du
end

"""Called after each full time step (not RK substage) for solution limiting."""
function post_step!(state, method::AbstractCapturingMethod, eq)
    return nothing
end

# Fine-grained sensor / dissipation API (AV family; optional composition)

function sense!(σ, sensor::AbstractShockSensor, u_work, state, eq)
    fill!(σ, zero(eltype(σ)))
    return σ
end

function sense!(σ, ::NullSensor, u_work, state, eq)
    fill!(σ, zero(eltype(σ)))
    return σ
end

function apply_dissipation!(du, dissip::AbstractDissipationOperator, σ, u_work, state, eq)
    return du
end

function apply_dissipation!(du, ::NullDissipation, σ, u_work, state, eq)
    return du
end

# ---------------------------------------------------------------------------
# Method registry (inventable names; further methods via src/methods/Registry.jl)
# ---------------------------------------------------------------------------

const METHOD_REGISTRY = Dict{String,Function}(
    "null" => (; kwargs...) -> NullCapturing(),
    # "persson_av" registered after PerssonAV.jl is loaded
)

function register_method!(name::AbstractString, factory::Function)
    METHOD_REGISTRY[String(name)] = factory
    return nothing
end

function get_capturing_method(name::AbstractString; kwargs...)
    key = String(name)
    haskey(METHOD_REGISTRY, key) ||
        error("Unknown capturing method \"$key\". Known: $(join(sort(collect(keys(METHOD_REGISTRY))), ", "))")
    return METHOD_REGISTRY[key](; kwargs...)
end
