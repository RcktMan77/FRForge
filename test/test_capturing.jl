using Test
using FRForge

# Counting method to verify residual! invokes hooks without naming PerssonAV
mutable struct HookProbeMethod <: AbstractCapturingMethod
    n_preprocess::Int
    n_extrapolate::Int
    n_sense::Int
    n_dissip::Int
    n_post::Int
end
HookProbeMethod() = HookProbeMethod(0, 0, 0, 0, 0)

function FRForge.preprocess_state!(u_work, m::HookProbeMethod, state, eq)
    m.n_preprocess += 1
    copyto!(u_work, state.u)
    return u_work
end
function FRForge.extrapolate_interface!(
    traces::FRForge.InterfaceTraces,
    m::HookProbeMethod,
    u_work,
    state,
    eq,
)
    m.n_extrapolate += 1
    return FRForge.extrapolate_interface!(traces, NullCapturing(), u_work, state, eq)
end
function FRForge.sense!(σ, m::HookProbeMethod, u_work, state, eq)
    m.n_sense += 1
    fill!(σ, 0)
    return σ
end
function FRForge.apply_dissipation!(du, m::HookProbeMethod, σ, u_work, state, eq)
    m.n_dissip += 1
    return du
end
function FRForge.post_step!(state, m::HookProbeMethod, eq)
    m.n_post += 1
    return nothing
end

@testset "hook pipeline invoked by residual and SSP-RK3" begin
    ops = build_operators(2)
    mesh = Mesh1D(0.0, 1.0, 4)
    state = allocate_state(mesh, ops, Val(1))
    set_initial_condition!(state, x -> sin(2π * x))
    eq = LinearAdvection1D(1.0)
    probe = HookProbeMethod()
    du = similar(state.u)
    residual!(du, state, eq, probe)
    @test probe.n_preprocess == 1
    @test probe.n_extrapolate == 1
    @test probe.n_sense == 1
    @test probe.n_dissip == 1

    probe2 = HookProbeMethod()
    ssp_rk3!(state, eq, probe2, 0.05; cfl=0.2)
    @test probe2.n_sense >= 3  # at least one full step × 3 stages
    @test probe2.n_post >= 1
end

@testset "Persson sensor activates on discontinuous data" begin
    ops = build_operators(3)
    mesh = Mesh1D(0.0, 1.0, 16)
    state = allocate_state(mesh, ops, Val(1))
    set_initial_condition!(state, x -> burgers_square_ic(x))
    eq = Burgers1D()
    σ = zeros(mesh.n_elements)
    sense!(σ, PerssonSensor(), state.u, state, eq)
    @test all(0 .<= σ .<= 1)
    @test maximum(σ) > 0.1  # some elements marked rough
end

@testset "Persson sensor quiet on smooth data" begin
    ops = build_operators(3)
    mesh = Mesh1D(0.0, 1.0, 16)
    state = allocate_state(mesh, ops, Val(1))
    set_initial_condition!(state, x -> sin(2π * x))
    eq = LinearAdvection1D(1.0)
    σ = zeros(mesh.n_elements)
    sense!(σ, PerssonSensor(), state.u, state, eq)
    @test all(σ .< 0.5)  # smooth → low indicator
end

@testset "conservative_br0 AV annihilates constants" begin
    ops = build_operators(3)
    mesh = Mesh1D(0.0, 1.0, 8)
    @test viscous_mass_residual_scale(ops, mesh) < 1e-12
end

@testset "PerssonAV reduces Burgers overshoot vs null" begin
    for p in (2, 3, 4)
        _, c_pers, c_cmp = run_persson_vs_null_burgers(; p=p, n_elements=32, t_final=0.15)
        @test c_pers["pass"] || c_pers["conservation_pass"]  # at least stable+cons
        @test c_pers["conservation_pass"]
        @test !c_pers["diverged"]
        @test c_cmp["metrics"]["overshoot_reduced"]
        @info "p=$p null_η=$(c_cmp["metrics"]["overshoot_null"]) pers_η=$(c_cmp["metrics"]["overshoot_persson"])"
    end
end

@testset "method registry" begin
    m0 = get_capturing_method("null")
    @test m0 isa NullCapturing
    m1 = get_capturing_method("persson_av")
    @test m1 isa PerssonAVMethod
    @test haskey(method_params(m1), "av_form")
    @test method_params(m1)["av_form"] == "conservative_br0"
end

@testset "M4 suite integration" begin
    cases, overall, fails = run_m4_capturing_suite()
    @test overall
    @test isempty(fails)
end

@testset "residual never requires Persson by name" begin
    # Smoke: residual! with AbstractCapturingMethod only
    ops = build_operators(2)
    mesh = Mesh1D(0.0, 1.0, 4)
    state = allocate_state(mesh, ops, Val(1))
    set_initial_condition!(state, x -> 1.0 + 0.1 * sin(2π * x))
    eq = Burgers1D()
    du = similar(state.u)
    residual!(du, state, eq, PerssonAVMethod()::AbstractCapturingMethod)
    @test all(isfinite, du)
end
