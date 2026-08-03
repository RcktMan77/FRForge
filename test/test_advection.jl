using Test
using FRForge

@testset "advection residual finite" begin
    p = 3
    ops = build_operators(p)
    mesh = Mesh1D(0.0, 1.0, 8; left_bc=PeriodicBC(), right_bc=PeriodicBC())
    state = allocate_state(mesh, ops, Val(1))
    set_initial_condition!(state, x -> sin(2π * x))
    eq = LinearAdvection1D(1.0)
    du = similar(state.u)
    residual!(du, state, eq, NullCapturing())
    @test all(isfinite, du)
    # Non-trivial residual for non-constant solution
    @test maximum(abs, du) > 1e-8
end

@testset "advection conservation p=2,3,4" begin
    for p in (2, 3, 4)
        c = run_advection_conservation(; p=p, n_elements=32, t_final=1.0, cfl=0.2)
        @test c["pass"]
        @test c["conservation_pass"]
        @test c["conservation_residual"] < 1e-10
        @test c["metrics"]["mass_abs_change"] < 1e-12
        @test !c["diverged"]
    end
end

@testset "advection smooth order p=2,3,4" begin
    for p in (2, 3, 4)
        nlist = p >= 4 ? [4, 8, 16] : [8, 16, 32, 64]
        c = run_advection_smooth_order(; p=p, n_elements_list=nlist, t_final=1.0, cfl=0.2)
        @test !c["diverged"]
        @test c["order_pass"]
        @test c["formal_order"] == p + 1
        # Each observed order within ~0.3 of formal
        for q in c["observed_orders"]
            @test isfinite(q)
            @test q >= (p + 1) - 0.3
        end
        @info "p=$p observed_orders=$(c["observed_orders"]) L2=$(c["l2_errors"])"
    end
end

@testset "M1 suite integration" begin
    cases, overall, fails = run_m1_advection_suite()
    @test overall
    @test isempty(fails)
    @test length(cases) >= 6  # 3 order + 3 conservation
end
