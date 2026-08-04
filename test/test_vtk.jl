using Test
using FRForge

@testset "VTK Lagrange line node order" begin
    # p=1: only endpoints
    n1 = vtk_lagrange_line_nodes(1)
    @test n1 ≈ [-1.0, 1.0]
    # p=2: left, right, mid
    n2 = vtk_lagrange_line_nodes(2)
    @test n2[1] ≈ -1
    @test n2[2] ≈ 1
    @test n2[3] ≈ 0
    # p=3: left, right, then two interiors
    n3 = vtk_lagrange_line_nodes(3)
    @test length(n3) == 4
    @test n3[1] ≈ -1
    @test n3[2] ≈ 1
    @test n3[3] ≈ -1 + 2 / 3
    @test n3[4] ≈ -1 + 4 / 3
end

@testset "point count oracle" begin
    np, nc, cl = vtk_point_counts_1d(5, 3)
    @test np == 5 * 4
    @test nc == 5
    @test cl == 4
end

@testset "GL to equi interpolation reproduces polynomials" begin
    for p in 1:4
        ops = build_operators(p)
        I = gl_to_equi_interp(ops)
        ξe = vtk_lagrange_line_nodes(p)
        # u = ξ^k exact for k ≤ p
        for k in 0:p
            u_gl = ops.ξ .^ k
            u_eq = I * u_gl
            @test u_eq ≈ ξe .^ k atol = 1e-12
        end
    end
end

@testset "write VTU advection p=2" begin
    mktempdir() do dir
        p = 2
        ops = build_operators(p)
        mesh = Mesh1D(0.0, 1.0, 4; left_bc=PeriodicBC(), right_bc=PeriodicBC())
        state = allocate_state(mesh, ops, Val(1))
        set_initial_condition!(state, x -> sin(2π * x))
        eq = LinearAdvection1D(1.0)
        path = joinpath(dir, "adv.vtu")
        write_vtu_high_order(path, state, eq)
        @test isfile(path)
        info = parse_vtu_basic(path)
        np, nc, _ = vtk_point_counts_1d(4, p)
        @test info.n_points == np
        @test info.n_cells == nc
        @test all(t -> t == VTK_LAGRANGE_LINE, info.types)
        @test length(info.types) == nc
        @test length(info.connectivity) == np
        # HO content: field not constant linear chord
        @test occursin("Name=\"u\"", info.text)
        # p=2 connectivity first cell: 0 1 2
        @test info.connectivity[1:3] == [0, 1, 2]
    end
end

@testset "write VTU Euler density wave" begin
    mktempdir() do dir
        p = 2
        eq = Euler1D(1.4)
        ops = build_operators(p)
        mesh = Mesh1D(0.0, 1.0, 3; left_bc=PeriodicBC(), right_bc=PeriodicBC())
        state = allocate_state(mesh, ops, Val(3))
        set_initial_condition!(state, x -> euler_density_wave_conserved(eq, x, 0.0))
        path = joinpath(dir, "euler.vtu")
        write_vtu_high_order(path, state, eq)
        info = parse_vtu_basic(path)
        @test info.n_points == 3 * 3
        @test all(==(VTK_LAGRANGE_LINE), info.types)
        txt = info.text
        for name in ("rho", "u", "p", "rho_u", "E")
            @test occursin("Name=\"$name\"", txt)
        end
    end
end

@testset "discontinuous points allow jumps" begin
    # Two elements, discontinuous IC — duplicated nodes at interface
    p = 1
    ops = build_operators(p)
    mesh = Mesh1D(0.0, 1.0, 2)
    state = allocate_state(mesh, ops, Val(1))
    set_initial_condition!(state, x -> x < 0.5 ? 1.0 : 0.0)
    eq = Burgers1D()
    mktempdir() do dir
        path = joinpath(dir, "jump.vtu")
        write_vtu_high_order(path, state, eq)
        info = parse_vtu_basic(path)
        @test info.n_points == 4  # 2 elements × 2 nodes (duplicated interface)
        @test info.n_cells == 2
    end
end

@testset "2D VTU capturing diagnostics (docs path)" begin
    # Optional sensor/av fields — documentation only, not invent default
    c, state, eq = run_euler2d_riemann(;
        p=1,
        nx=6,
        ny=6,
        t_final=0.02,
        cfl=0.05,
        method=PerssonAVMethod(; c_av=0.1),
        method_name="persson_av",
    )
    @test !c["diverged"]
    method = PerssonAVMethod(; c_av=0.1)
    σ, ε = compute_capturing_diagnostics_2d(state, eq, method)
    @test length(σ) == state.mesh.n_elements
    @test all(σ .>= 0)
    mktempdir() do dir
        path = joinpath(dir, "riemann_diag.vtu")
        write_vtu_high_order_with_capturing(path, state, eq, method)
        txt = read(path, String)
        @test occursin("Name=\"rho\"", txt)
        @test occursin("Name=\"sensor\"", txt)
        @test occursin("Name=\"av\"", txt)
        # default writer still omits diagnostics
        path2 = joinpath(dir, "plain.vtu")
        write_vtu_high_order(path2, state, eq)
        txt2 = read(path2, String)
        @test occursin("Name=\"rho\"", txt2)
        @test !occursin("Name=\"sensor\"", txt2)
    end
end
