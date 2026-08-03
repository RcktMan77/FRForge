# SSP-RK2 (Heun / explicit trapezoidal) — second-order strong-stability-preserving RK.
#
# Stages (Gottlieb–Shu):
#   u⁽¹⁾ = uⁿ + Δt L(uⁿ)
#   uⁿ⁺¹ = ½ uⁿ + ½ u⁽¹⁾ + ½ Δt L(u⁽¹⁾)
#
# CFL: use the same form as SSP-RK3 (default CFL ≈ 0.2). Prefer CFL ≤ that used
# for SSP-RK3; for high-order spatial studies prefer fixed small Δt because the
# scheme is only second-order in time.

"""
    ssp_rk2_step!(state, eq, method, Δt; du, u0, u1) -> Symbol
"""
function ssp_rk2_step!(
    state::SolutionState{T,Neq},
    eq::AbstractEquation{Neq},
    method::AbstractCapturingMethod,
    Δt::T;
    du=similar(state.u),
    u0=similar(state.u),
    u1=similar(state.u),
) where {T,Neq}
    copyto!(u0, state.u)
    t0 = state.t

    residual!(du, state, eq, method)
    has_nonfinite(du) && return :diverged
    @. state.u = u0 + Δt * du
    state.t = t0 + Δt

    residual!(du, state, eq, method)
    has_nonfinite(du) && return :diverged
    @. u1 = state.u
    @. state.u = (T(1) / T(2)) * u0 + (T(1) / T(2)) * u1 + (T(1) / T(2)) * Δt * du
    state.t = t0 + Δt

    has_nonfinite(state.u) && return :diverged
    post_step!(state, method, eq)
    return :ok
end

function ssp_rk2_step!(
    state::SolutionState2D{T,Neq},
    eq::AbstractEquation{Neq},
    method::AbstractCapturingMethod,
    Δt::T;
    du=similar(state.u),
    u0=similar(state.u),
    u1=similar(state.u),
) where {T,Neq}
    copyto!(u0, state.u)
    t0 = state.t
    residual!(du, state, eq, method)
    has_nonfinite(du) && return :diverged
    @. state.u = u0 + Δt * du
    state.t = t0 + Δt
    residual!(du, state, eq, method)
    has_nonfinite(du) && return :diverged
    @. u1 = state.u
    @. state.u = (T(1) / T(2)) * u0 + (T(1) / T(2)) * u1 + (T(1) / T(2)) * Δt * du
    state.t = t0 + Δt
    has_nonfinite(state.u) && return :diverged
    post_step!(state, method, eq)
    return :ok
end

"""
    ssp_rk2!(state, eq, method, t_final; cfl=0.2, dt=nothing, max_steps=10^7)
"""
function ssp_rk2!(
    state::SolutionState{T,Neq},
    eq::AbstractEquation{Neq},
    method::AbstractCapturingMethod,
    t_final::Real;
    cfl::Real=0.2,
    dt::Union{Nothing,Real}=nothing,
    max_steps::Int=10^7,
) where {T,Neq}
    t_final_T = T(t_final)
    du = similar(state.u)
    u0 = similar(state.u)
    u1 = similar(state.u)
    fixed_dt = dt === nothing ? nothing : T(dt)
    n_steps = 0
    while state.t < t_final_T - 10 * eps(T)
        step_dt = fixed_dt === nothing ? compute_dt(state, eq; cfl=cfl) : fixed_dt
        if state.t + step_dt > t_final_T
            step_dt = t_final_T - state.t
        end
        status = ssp_rk2_step!(state, eq, method, step_dt; du=du, u0=u0, u1=u1)
        n_steps += 1
        if status != :ok
            return (status=status, n_steps=n_steps, t=state.t)
        end
        if n_steps >= max_steps
            return (status=:max_steps, n_steps=n_steps, t=state.t)
        end
    end
    return (status=:ok, n_steps=n_steps, t=state.t)
end

function ssp_rk2!(
    state::SolutionState2D{T,Neq},
    eq::AbstractEquation{Neq},
    method::AbstractCapturingMethod,
    t_final::Real;
    cfl::Real=0.2,
    dt::Union{Nothing,Real}=nothing,
    max_steps::Int=10^7,
) where {T,Neq}
    t_final_T = T(t_final)
    du = similar(state.u)
    u0 = similar(state.u)
    u1 = similar(state.u)
    fixed_dt = dt === nothing ? nothing : T(dt)
    n_steps = 0
    while state.t < t_final_T - 10 * eps(T)
        step_dt = fixed_dt === nothing ? compute_dt(state, eq; cfl=cfl) : fixed_dt
        if state.t + step_dt > t_final_T
            step_dt = t_final_T - state.t
        end
        status = ssp_rk2_step!(state, eq, method, step_dt; du=du, u0=u0, u1=u1)
        n_steps += 1
        if status != :ok
            return (status=status, n_steps=n_steps, t=state.t)
        end
        if n_steps >= max_steps
            return (status=:max_steps, n_steps=n_steps, t=state.t)
        end
    end
    return (status=:ok, n_steps=n_steps, t=state.t)
end

function ssp_rk2!(state, eq, t_final; kwargs...)
    return ssp_rk2!(state, eq, NullCapturing(), t_final; kwargs...)
end
