using Test
using FRForge

@testset "Burgers residual finite" begin
    p = 3
    ops = build_operators(p)
    mesh = Mesh1D(0.0, 1.0, 16; left_bc = PeriodicBC(), right_bc = PeriodicBC())
    state = allocate_state(mesh, ops, Val(1))
    set_initial_condition!(state, x -> burgers_square_ic(x))
    eq = Burgers1D()
    du = similar(state.u)
    residual!(du, state, eq, NullCapturing())
    @test all(isfinite, du)
    @test maximum(abs, du) > 1e-8
end

@testset "Rusanov flux for Burgers" begin
    eq = Burgers1D()
    # Equal states → physical flux
    f = numerical_flux(eq, [1.0], [1.0])
    @test f[1] ≈ 0.5 atol = 1e-14
    # Upwind-ish dissipation present when states differ
    f2 = numerical_flux(eq, [1.0], [0.0])
    @test f2[1] ≈ 0.5 * (0.5 + 0.0) - 0.5 * 1.0 * (0.0 - 1.0)  # 0.25 + 0.5 = 0.75
end

@testset "burgers conservation p=2,3,4" begin
    for p in (2, 3, 4)
        c = run_burgers_conservation(; p = p, n_elements = 32, t_final = 0.15, cfl = 0.2)
        @test c["pass"]
        @test c["conservation_pass"]
        @test c["conservation_residual"] < 1e-10
        @test !c["diverged"]
    end
end

@testset "burgers HO oscillations p=2,3,4" begin
    for p in (2, 3, 4)
        c = run_burgers_oscillation(; p = p, n_elements = 32, t_final = 0.15, cfl = 0.2)
        @test !c["diverged"]
        @test c["conservation_pass"]
        @test c["metrics"]["oscillations_present"]
        @test c["overshoot"] >= 0.02
        @test c["pass"]  # success means oscillations were demonstrated
        @info "p=$p overshoot=$(c["overshoot"]) u∈[$(c["metrics"]["u_min"]), $(c["metrics"]["u_max"])]"
    end
end
