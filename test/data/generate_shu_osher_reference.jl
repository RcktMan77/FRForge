#!/usr/bin/env julia
# Generate test/data/shu_osher_ref.csv and update README hash.
# Usage: julia --project=. test/data/generate_shu_osher_reference.jl

using FRForge
using SHA

const OUT = joinpath(@__DIR__, "shu_osher_ref.csv")
const README = joinpath(@__DIR__, "README.md")

function generate(; p=1, ne=200, c_av=0.5, cfl=0.1, t_final=1.8)
    eq = Euler1D(1.4)
    ops = build_operators(p)
    mesh = Mesh1D(0.0, 10.0, ne; left_bc=TransmissiveBC(), right_bc=TransmissiveBC())
    state = allocate_state(mesh, ops, Val(3))
    set_initial_condition!(state, x -> FRForge.shu_osher_ic(eq, x))
    method = PerssonAVMethod(c_av=c_av)
    result = ssp_rk3!(state, eq, method, t_final; cfl=cfl)
    result.status == :ok || error("generation failed: $(result.status)")
    positivity_ok(eq, state) || error("generation lost positivity")
    x, ρ = sample_solution_1d(state, eq; component=:density)
    open(OUT, "w") do io
        println(io, "# Shu-Osher density reference")
        println(io, "# PerssonAV p=$p Ne=$ne c_av=$c_av t=$t_final domain=[0,10] CFL=$cfl")
        println(io, "# columns: x, rho")
        for i in eachindex(x)
            println(io, x[i], ",", ρ[i])
        end
    end
    h = bytes2hex(sha256(read(OUT)))
    open(README, "w") do io
        write(
            io,
            """
# Verification reference data

## shu_osher_ref.csv

| Field | Value |
|-------|-------|
| Problem | Shu–Osher on [0,10], interface x=1, t=$t_final |
| Generator | PerssonAV FR (p=$p, Ne=$ne, c_av=$c_av, CFL=$cfl, SSP-RK3) |
| Columns | `x, rho` |
| SHA-256 | `$h` |
| ρ range | [$(minimum(ρ)), $(maximum(ρ))] |

### Policy (design OQ1)

1. Self-generated fine-grid capturing run (not external-only).
2. Self-convergence: compare Ne=100/200/400 shock locus and density extrema.
3. Sanity: primary shock near x≈2.4; high-frequency post-shock waves; ρ ~ O(1)–O(4).
4. NullCapturing is often unstable; frozen reference uses trusted Persson AV.
5. Regenerate: `julia --project=. test/data/generate_shu_osher_reference.jl`
""",
        )
    end
    println("Wrote $OUT")
    println("SHA-256 $h")
    println("ρ ∈ [$(minimum(ρ)), $(maximum(ρ))]  n=$(length(x))")
    return h
end

# Optional coarse self-convergence diagnostics
function self_convergence_check()
    println("Self-convergence (p=1, c_av=0.5):")
    for ne in (100, 200)
        eq = Euler1D(1.4)
        ops = build_operators(1)
        mesh = Mesh1D(0.0, 10.0, ne; left_bc=TransmissiveBC(), right_bc=TransmissiveBC())
        state = allocate_state(mesh, ops, Val(3))
        set_initial_condition!(state, x -> FRForge.shu_osher_ic(eq, x))
        ssp_rk3!(state, eq, PerssonAVMethod(c_av=0.5), 1.8; cfl=0.1)
        x, ρ = sample_solution_1d(state, eq; component=:density)
        # shock locus ~ max |dρ|
        imax = argmax(abs.(diff(ρ)) ./ max.(abs.(diff(x)), 1e-30))
        println("  Ne=$ne shock_x≈$(x[imax]) ρ∈[$(minimum(ρ)), $(maximum(ρ))]")
    end
end

generate()
self_convergence_check()
