using Test
using FRForge
using LinearAlgebra

@testset "SchemeConfig defaults frozen" begin
    @test DEFAULT_SCHEME.points === :gl
    @test DEFAULT_SCHEME.flux === :rusanov
    @test DEFAULT_SCHEME.time === :ssp_rk3
    d = scheme_dict(DEFAULT_SCHEME)
    @test d["points"] == "GL"
    @test d["flux"] == "Rusanov"
    @test d["time"] == "SSP-RK3"
    s = parse_scheme(; points="gll", flux="hllc", time="ssp-rk2")
    @test s.points === :gll
    @test s.flux === :hllc
    @test s.time === :ssp_rk2
end

@testset "GLL nodes/weights and operators" begin
    for n in 2:6
        ξ, w = gauss_lobatto_legendre_nodes_weights(n)
        @test length(ξ) == n
        @test ξ[1] ≈ -1
        @test ξ[end] ≈ 1
        @test sum(w) ≈ 2 atol = 1e-12
        @test issorted(ξ)
    end
    for p in 1:4
        ops = build_operators(p; points=:gll)
        @test ops.points === :gll
        @test n_points(ops) == p + 1
        @test ops.ξ[1] ≈ -1
        @test ops.ξ[end] ≈ 1
        # poly differentiation exactness
        for k in 0:p
            u = ops.ξ .^ k
            du = ops.D * u
            duex = k == 0 ? zero.(ops.ξ) : (k .* ops.ξ .^ (k - 1))
            @test du ≈ duex atol = 1e-10 rtol = 1e-10
        end
        @test sum(ops.ℓ_L) ≈ 1 atol = 1e-12
        @test sum(ops.ℓ_R) ≈ 1 atol = 1e-12
        # endpoints are nodes → unit Lagrange basis
        @test ops.ℓ_L[1] ≈ 1 atol = 1e-12
        @test ops.ℓ_R[end] ≈ 1 atol = 1e-12
    end
    # GL remains default
    ops_gl = build_operators(3)
    @test ops_gl.points === :gl
    @test abs(ops_gl.ξ[1] + 1) > 1e-8  # interior nodes
end

@testset "GLL advection order and conservation" begin
    p = 3
    ops = build_operators(p; points=:gll)
    scheme = SchemeConfig(; points=:gll)
    eq = LinearAdvection1D(1.0)
    mesh = Mesh1D(0.0, 1.0, 16; left_bc=PeriodicBC(), right_bc=PeriodicBC())
    state = allocate_state(mesh, ops, Val(1); scheme=scheme)
    @test state.scheme.points === :gll
    set_initial_condition!(state, x -> sin(2π * x))
    M0 = discrete_mass(state, 1)
    # fixed small dt for spatial order (SSP-RK3 default time)
    result = integrate!(state, eq, NullCapturing(), 1.0; dt=1e-4)
    @test result.status == :ok
    @test abs(discrete_mass(state, 1) - M0) < 1e-10
    err = l2_error(state, x -> sin(2π * (x - 1.0)), 1)
    @test err < 1e-3

    # modest order study
    errs = Float64[]
    for Ne in (8, 16, 32)
        mesh = Mesh1D(0.0, 1.0, Ne; left_bc=PeriodicBC(), right_bc=PeriodicBC())
        st = allocate_state(mesh, ops, Val(1); scheme=scheme)
        set_initial_condition!(st, x -> sin(2π * x))
        dt = 0.05 / (Ne * (2p + 1))
        r = integrate!(st, eq, NullCapturing(), 1.0; dt=dt)
        @test r.status == :ok
        push!(errs, l2_error(st, x -> sin(2π * (x - 1.0)), 1))
    end
    # observed order between 8→16 and 16→32
    o1 = log(errs[1] / errs[2]) / log(2)
    o2 = log(errs[2] / errs[3]) / log(2)
    @test o1 > p - 0.5
    @test o2 > p - 0.5
end

@testset "HLLC Sod interface flux" begin
    eq = Euler1D(1.4)
    uL = primitives_to_conserved(eq, 1.0, 0.0, 1.0)
    uR = primitives_to_conserved(eq, 0.125, 0.0, 0.1)
    f_r = rusanov_flux(eq, uL, uR)
    f_h = hllc_flux(eq, uL, uR)
    @test length(f_h) == 3
    @test all(isfinite, f_h)
    @test all(isfinite, f_r)
    # Stationary contact: equal pressure and velocity, density jump
    uLc = primitives_to_conserved(eq, 1.0, 0.0, 1.0)
    uRc = primitives_to_conserved(eq, 0.1, 0.0, 1.0)
    fh = hllc_flux(eq, uLc, uRc)
    # momentum flux = p for stationary fluid
    @test fh[2] ≈ 1.0 atol = 1e-10
end

@testset "HLLC Euler density wave (NullCapturing)" begin
    p = 2
    eq = Euler1D(1.4)
    ops = build_operators(p; points=:gl)
    scheme = SchemeConfig(; flux=:hllc)
    mesh = Mesh1D(0.0, 1.0, 32; left_bc=PeriodicBC(), right_bc=PeriodicBC())
    state = allocate_state(mesh, ops, Val(3); scheme=scheme)
    set_initial_condition!(state, x -> euler_density_wave_conserved(eq, x, 0.0))
    M0 = discrete_mass(state, 1)
    result = integrate!(state, eq, NullCapturing(), 0.1; cfl=0.1)
    @test result.status == :ok
    @test positivity_ok(eq, state.u)
    @test abs(discrete_mass(state, 1) - M0) < 1e-10
    err = l2_error(state, x -> euler_density_wave_conserved(eq, x, 0.1), 1)
    @test err < 5e-3
end

@testset "SSP-RK2 advection smoke" begin
    p = 2
    ops = build_operators(p)
    scheme = SchemeConfig(; time=:ssp_rk2)
    eq = LinearAdvection1D(1.0)
    mesh = Mesh1D(0.0, 1.0, 16; left_bc=PeriodicBC(), right_bc=PeriodicBC())
    state = allocate_state(mesh, ops, Val(1); scheme=scheme)
    set_initial_condition!(state, x -> sin(2π * x))
    M0 = discrete_mass(state, 1)
    result = integrate!(state, eq, NullCapturing(), 1.0; dt=5e-4)
    @test result.status == :ok
    @test abs(discrete_mass(state, 1) - M0) < 1e-10
    err = l2_error(state, x -> sin(2π * (x - 1.0)), 1)
    @test err < 5e-3
    @test occursin("SSP-RK2", time_cfl_guidance(:ssp_rk2))
end

@testset "scheme combo GLL×HLLC×SSP-RK2 smoke" begin
    p = 2
    scheme = SchemeConfig(; points=:gll, flux=:hllc, time=:ssp_rk2)
    ops = build_operators(p; points=scheme.points)
    eq = Euler1D(1.4)
    mesh = Mesh1D(0.0, 1.0, 24; left_bc=PeriodicBC(), right_bc=PeriodicBC())
    state = allocate_state(mesh, ops, Val(3); scheme=scheme)
    set_initial_condition!(state, x -> euler_density_wave_conserved(eq, x, 0.0))
    result = integrate!(state, eq, NullCapturing(), 0.05; cfl=0.1)
    @test result.status == :ok
    @test positivity_ok(eq, state.u)
end

@testset "report includes scheme" begin
    mktempdir() do dir
        path = joinpath(dir, "r.json")
        r = write_report_skeleton(
            path;
            command="test",
            suite="smoke",
            scheme=SchemeConfig(; points=:gll, flux=:hllc, time=:ssp_rk2),
        )
        @test r["scheme"]["points"] == "GLL"
        @test r["scheme"]["flux"] == "HLLC"
        @test r["scheme"]["time"] == "SSP-RK2"
        @test isempty(validate_report_keys(r))
    end
end
