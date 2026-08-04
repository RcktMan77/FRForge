# Scheme-aware time integration dispatcher.

"""
    integrate!(state, eq, method, t_final; cfl, dt, max_steps, time=nothing)

Advance `state` to `t_final` using `state.scheme.time` (or override `time`).
Defaults to SSP-RK3 (frozen invent scheme).
"""
function integrate!(
    state,
    eq,
    method::AbstractCapturingMethod,
    t_final::Real;
    cfl::Real=0.2,
    dt::Union{Nothing,Real}=nothing,
    max_steps::Int=10^7,
    time::Union{Nothing,Symbol}=nothing,
    progress_every::Int=0,
    progress_label::AbstractString="",
)
    tsym = something(time, state.scheme.time)
    if tsym === :ssp_rk2
        return ssp_rk2!(state, eq, method, t_final; cfl=cfl, dt=dt, max_steps=max_steps)
    elseif tsym === :ssp_rk3
        return ssp_rk3!(
            state,
            eq,
            method,
            t_final;
            cfl=cfl,
            dt=dt,
            max_steps=max_steps,
            progress_every=progress_every,
            progress_label=progress_label,
        )
    else
        error("Unknown time scheme $tsym (use :ssp_rk3 or :ssp_rk2)")
    end
end

function integrate!(state, eq, t_final; kwargs...)
    return integrate!(state, eq, NullCapturing(), t_final; kwargs...)
end
