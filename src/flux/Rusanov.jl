# Rusanov (local Lax–Friedrichs) numerical flux.

"""
    rusanov_flux(eq, uL, uR) -> numerical flux

Local Lax–Friedrichs:
  f̂ = ½(F(u⁻)+F(u⁺)) − ½ λ_max (u⁺ − u⁻)
with λ_max from `max_wave_speed(eq, uL, uR)`.
"""
function rusanov_flux(eq::AbstractEquation{Neq}, uL::AbstractVector{T}, uR::AbstractVector{T}) where {T,Neq}
    fL = physical_flux(eq, uL)
    fR = physical_flux(eq, uR)
    λ = max_wave_speed(eq, uL, uR)
    out = Vector{T}(undef, Neq)
    @inbounds for c in 1:Neq
        fl = fL isa AbstractVector ? fL[c] : T(fL)
        fr = fR isa AbstractVector ? fR[c] : T(fR)
        out[c] = T(0.5) * (fl + fr) - T(0.5) * T(λ) * (uR[c] - uL[c])
    end
    return out
end

function rusanov_flux(eq::AbstractEquation{1}, uL::Number, uR::Number)
    T = typeof(float(uL))
    fL = physical_flux(eq, uL)
    fR = physical_flux(eq, uR)
    λ = max_wave_speed(eq, uL, uR)
    return T(0.5) * (T(fL) + T(fR)) - T(0.5) * T(λ) * (T(uR) - T(uL))
end
