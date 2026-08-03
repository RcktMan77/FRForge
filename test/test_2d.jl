using Test
using FRForge

@testset "Mesh2D basics" begin
    mesh = Mesh2D(0.0, 1.0, 0.0, 1.0, 4, 3)
    @test mesh.n_elements == 12
    @test element_index(mesh, 2, 1) == 2
    @test element_coords(mesh, 5) == (1, 2)
end

@testset "2D advection residual finite" begin
    eq = LinearAdvection2D(1.0, 0.5)
    ops = build_operators(2)
    mesh = Mesh2D(0.0, 1.0, 0.0, 1.0, 4, 4)
    state = allocate_state(mesh, ops, Val(1))
    set_initial_condition!(state, (x, y) -> sin(2π * x) * sin(2π * y))
    du = similar(state.u)
    residual!(du, state, eq, NullCapturing())
    @test all(isfinite, du)
    @test maximum(abs, du) > 1e-8
end

@testset "2D advection order p=2" begin
    c = run_advection2d_smooth_order(; p=2, n_list=[4, 8, 16], t_final=1.0, cfl=0.12)
    @test !c["diverged"]
    @test c["order_pass"]
    @info "2D adv orders=$(c["observed_orders"]) L2=$(c["l2_errors"])"
end

@testset "2D Euler order p=2" begin
    c = run_euler2d_smooth_order(; p=2, n_list=[16, 32], t_final=0.25, cfl=0.06)
    @test !c["diverged"]
    @test c["order_pass"]
    @test c["positivity_ok"]
    @info "2D Euler orders=$(c["observed_orders"])"
end

@testset "2D discontinuous Euler runs" begin
    c, state, eq = run_euler2d_discontinuous(; p=1, nx=12, ny=12, t_final=0.05)
    @test c["pass"]
    @test c["positivity_ok"]
end

@testset "2D Persson AV sensor and viscous mass" begin
    ops = build_operators(2)
    mesh = Mesh2D(0.0, 1.0, 0.0, 1.0, 4, 4)
    @test viscous_mass_residual_scale_2d(ops, mesh) < 1e-10

    eq = Euler2D(1.4)
    state = allocate_state(mesh, ops, Val(4))
    # Interior-element jump (not only face-aligned) so modal sensor activates
    set_initial_condition!(
        state,
        (x, y) ->
            (x + y < 0.7) ? primitives_to_conserved(eq, 1.0, 0.0, 0.0, 1.0) :
            primitives_to_conserved(eq, 0.125, 0.0, 0.0, 0.1),
    )
    method = PerssonAVMethod(; c_av=0.1)
    σ = zeros(mesh.n_elements)
    sense!(σ, method, state.u, state, eq)
    @test all(0 .<= σ .<= 1)
    @test maximum(σ) > 0.01  # under-resolved jump should trigger somewhere
    du = zeros(size(state.u)...)
    residual!(du, state, eq, method)
    @test all(isfinite, du)
end

@testset "2D Persson AV discontinuous CI-light" begin
    c, _, _ = run_euler2d_discontinuous(;
        p=1,
        nx=12,
        ny=12,
        t_final=0.05,
        cfl=0.08,
        method=PerssonAVMethod(; c_av=0.1),
        method_name="persson_av",
    )
    @test c["pass"]
    @test c["positivity_ok"]
    @test !c["diverged"]
end

@testset "2D Persson AV smooth stability" begin
    c = run_euler2d_smooth_order(;
        p=2,
        n_list=[8, 16],
        t_final=0.15,
        cfl=0.05,
        method=PerssonAVMethod(; c_av=0.1),
        method_name="persson_av",
    )
    @test !c["diverged"]
    @test c["positivity_ok"]
    @test c["conservation_pass"]
end

@testset "P3.1 2D capturing suite" begin
    cases, overall, fails = run_p31_2d_capturing_suite()
    @test overall
    @test isempty(fails)
    @test length(cases) >= 3
end

@testset "2D VTU Lagrange quad" begin
    mktempdir() do dir
        eq = LinearAdvection2D(1.0, 0.0)
        ops = build_operators(2)
        mesh = Mesh2D(0.0, 1.0, 0.0, 1.0, 2, 2)
        state = allocate_state(mesh, ops, Val(1))
        set_initial_condition!(state, (x, y) -> x + y)
        path = joinpath(dir, "adv2d.vtu")
        write_vtu_high_order(path, state, eq)
        info = parse_vtu_basic(path)
        np, nc, cl = vtk_point_counts_2d(4, 2)
        @test info.n_points == np
        @test info.n_cells == nc
        @test all(==(VTK_LAGRANGE_QUAD), info.types)
        @test occursin("Name=\"u\"", info.text)
    end
end

@testset "M8 suite" begin
    cases, overall, fails, _, _ = run_m8_2d_suite()
    @test overall
    @test isempty(fails)
    @test length(cases) >= 3
end
