using Test
using FRForge

@testset "Euler primitives and flux" begin
    eq = Euler1D(1.4)
    U = primitives_to_conserved(eq, 1.0, 0.5, 1.0)
    ρ, u, p = conserved_to_primitives(eq, U)
    @test ρ ≈ 1.0
    @test u ≈ 0.5
    @test p ≈ 1.0
    f = physical_flux(eq, U)
    @test f[1] ≈ U[2]
    @test f[2] ≈ U[2] * u + p
    @test positivity_ok_state(eq, U)
end

@testset "Euler residual finite" begin
    eq = Euler1D(1.4)
    ops = build_operators(3)
    mesh = Mesh1D(0.0, 1.0, 8; left_bc=PeriodicBC(), right_bc=PeriodicBC())
    state = allocate_state(mesh, ops, Val(3))
    set_initial_condition!(state, x -> euler_density_wave_conserved(eq, x, 0.0))
    du = similar(state.u)
    residual!(du, state, eq, NullCapturing())
    @test all(isfinite, du)
    @test maximum(abs, du) > 1e-10
end

@testset "euler conservation p=2,3,4" begin
    for p in (2, 3, 4)
        c = run_euler_conservation(; p=p, n_elements=32, t_final=1.0)
        @test c["pass"]
        @test c["conservation_pass"]
        @test c["positivity_ok"]
        @test c["conservation_residual"] < 1e-10
    end
end

@testset "euler smooth order p=2,3,4" begin
    for p in (2, 3, 4)
        nlist = p == 2 ? [16, 32, 64] : (p == 3 ? [8, 16, 32] : [16, 32, 64])
        c = run_euler_smooth_order(; p=p, n_elements_list=nlist, t_final=1.0, cfl=0.1)
        @test !c["diverged"]
        @test c["order_pass"]
        @test c["formal_order"] == p + 1
        @test c["positivity_ok"]
        for q in c["observed_orders"]
            @test isfinite(q)
            @test q >= (p + 1) - 0.3
        end
        @info "euler p=$p orders=$(c["observed_orders"]) L2=$(c["l2_errors"])"
    end
end

@testset "BC freestream tests" begin
    bt = run_bc_transmissive_test()
    @test bt["pass"]
    @test bt["metrics"]["max_deviation_from_uniform"] < 1e-12

    bd = run_bc_dirichlet_test()
    @test bd["pass"]
    @test bd["metrics"]["max_deviation_from_uniform"] < 1e-12
end

