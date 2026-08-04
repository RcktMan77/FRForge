# Example inventable method: Persson AV with elevated c_av (demonstrates invent pathway).
# Agents should replace this pattern with structurally novel hooks under src/methods/.

"""
    ScaledPerssonMethod

Thin wrapper around `PerssonAVMethod` with a higher default artificial viscosity
coefficient. Exists so the invent loop can be exercised end-to-end against
`persson_av` without requiring a brand-new algorithm on day one.

For real research, implement new types that override
`preprocess_state!`, `extrapolate_interface!`, `sense!`, `apply_dissipation!`,
and/or `post_step!` instead of only retuning `c_av`.
"""
struct ScaledPerssonMethod{T} <: AbstractCapturingMethod
    inner::PerssonAVMethod{T}
end

function ScaledPerssonMethod(;
    κ::Real = 4.0,
    s0_factor::Real = -4.0,
    c_av::Real = 0.5,  # elevated vs baseline 0.1 so invent demo can clear δ=0.02
    av_form::AbstractString = "conservative_br0",
    T::Type = Float64,
)
    return ScaledPerssonMethod{T}(
        PerssonAVMethod(; κ = κ, s0_factor = s0_factor, c_av = c_av, av_form = av_form, T = T),
    )
end

# Forward all hooks to the inner Persson method (composition pattern)
function preprocess_state!(u_work, m::ScaledPerssonMethod, state, eq)
    return preprocess_state!(u_work, m.inner, state, eq)
end
function extrapolate_interface!(
    traces::InterfaceTraces,
    m::ScaledPerssonMethod,
    u_work,
    state,
    eq,
)
    return extrapolate_interface!(traces, m.inner, u_work, state, eq)
end
function numerical_flux_method(m::ScaledPerssonMethod, eq, uL, uR)
    return numerical_flux_method(m.inner, eq, uL, uR)
end
function sense!(σ, m::ScaledPerssonMethod, u_work, state, eq)
    return sense!(σ, m.inner, u_work, state, eq)
end
function apply_dissipation!(du, m::ScaledPerssonMethod, σ, u_work, state, eq)
    return apply_dissipation!(du, m.inner, σ, u_work, state, eq)
end
function post_step!(state, m::ScaledPerssonMethod, eq)
    return post_step!(state, m.inner, eq)
end

function method_params(m::ScaledPerssonMethod)
    p = method_params(m.inner)
    p["method"] = "scaled_persson"
    p["note"] = "example inventable method; elevated c_av over baseline persson_av"
    return p
end
