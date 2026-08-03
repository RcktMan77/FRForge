# SSP-RK3 (Shu–Osher) time integration.

"""
    compute_dt(state, eq; cfl=0.2) -> Δt

Δt = CFL * min_e (Δx_e / ((2p+1) λ_max)).
"""
function compute_dt(state::SolutionState{T}, eq; cfl::Real=0.2) where {T}
    mesh, ops = state.mesh, state.ops
    p = ops.p
    λ = max_wave_speed(eq, state.u)
    λ = max(λ, eps(T))
    dt = typemax(T)
    for e in 1:mesh.n_elements
        dt = min(dt, T(cfl) * mesh.Δx[e] / (T(2p + 1) * T(λ)))
    end
    return dt
end

"""
    ssp_rk3_step!(state, eq, method, Δt; du, u0, u1, u2)

One SSP-RK3 step (Shu–Osher):
  u¹ = uⁿ + Δt L(uⁿ)
  u² = 3/4 uⁿ + 1/4 u¹ + 1/4 Δt L(u¹)
  uⁿ⁺¹ = 1/3 uⁿ + 2/3 u² + 2/3 Δt L(u²)
then post_step!.
"""
function ssp_rk3_step!(
    state::SolutionState{T,Neq},
    eq::AbstractEquation{Neq},
    method::AbstractCapturingMethod,
    Δt::T;
    du=similar(state.u),
    u0=similar(state.u),
    u1=similar(state.u),
    u2=similar(state.u),
) where {T,Neq}
    copyto!(u0, state.u)
    t0 = state.t

    # Stage 1
    residual!(du, state, eq, method)
    has_nonfinite(du) && return :diverged
    @. state.u = u0 + Δt * du
    state.t = t0 + Δt

    # Stage 2
    residual!(du, state, eq, method)
    has_nonfinite(du) && return :diverged
    @. u1 = state.u
    @. state.u = (T(3) / T(4)) * u0 + (T(1) / T(4)) * u1 + (T(1) / T(4)) * Δt * du
    state.t = t0 + Δt / T(2)

    # Stage 3
    residual!(du, state, eq, method)
    has_nonfinite(du) && return :diverged
    @. u2 = state.u
    @. state.u = (T(1) / T(3)) * u0 + (T(2) / T(3)) * u2 + (T(2) / T(3)) * Δt * du
    state.t = t0 + Δt

    has_nonfinite(state.u) && return :diverged
    post_step!(state, method, eq)
    return :ok
end

"""
    ssp_rk3!(state, eq, method, t_final; cfl=0.2, dt=nothing, max_steps=10^7) -> NamedTuple

Integrate until `t_final` with SSP-RK3. Returns status and step count.

If `dt` is provided, use that fixed step (last step may be shorter).
Otherwise `dt` is recomputed each step from `cfl`.
"""
function ssp_rk3!(
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
    u2 = similar(state.u)
    fixed_dt = dt === nothing ? nothing : T(dt)

    n_steps = 0
    while state.t < t_final_T - 10 * eps(T)
        step_dt = fixed_dt === nothing ? compute_dt(state, eq; cfl=cfl) : fixed_dt
        if state.t + step_dt > t_final_T
            step_dt = t_final_T - state.t
        end
        status = ssp_rk3_step!(state, eq, method, step_dt; du=du, u0=u0, u1=u1, u2=u2)
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

function ssp_rk3!(state, eq, t_final; kwargs...)
    return ssp_rk3!(state, eq, NullCapturing(), t_final; kwargs...)
end
