using Test
using FRForge

@testset "Sod exact solver sanity" begin
    prob = SodProblem()
    ρL, uL, pL = sod_exact(prob, 0.2, 0.2)
    @test ρL ≈ 1.0 atol = 1e-10
    @test abs(uL) < 1e-10
    # Right of shock should be near right state
    ρR, uR, pR = sod_exact(prob, 0.95, 0.2)
    @test ρR ≈ 0.125 atol = 1e-3
end

@testset "Sod run metrics" begin
    c = run_sod(; p=2, n_elements=48, method=PerssonAVMethod(), method_name="persson_av")
    @test c["pass"]
    @test c["positivity_ok"]
    @test c["case_type"] == "discontinuous"
    @test c["excess_dissipation"] !== nothing
    @test isfinite(c["shock_thickness"])
    @test c["shock_thickness"] < 40  # localized, not whole domain
    @test c["overshoot"] < 0.5
    @test c["l1_error_vs_reference"] < 0.05
end

@testset "Shu-Osher run" begin
    c = run_shu_osher(;
        p=1,
        n_elements=80,
        cfl=0.1,
        method=PerssonAVMethod(c_av=0.4),
        method_name="persson_av",
    )
    @test c["pass"]
    @test c["positivity_ok"]
    @test c["case_type"] == "discontinuous"
    @test isfinite(c["overshoot"])
end

@testset "Shu-Osher reference file" begin
    path = joinpath(@__DIR__, "data", "shu_osher_ref.csv")
    @test isfile(path)
    x, ρ = FRForge.load_shu_osher_reference(path)
    @test length(x) == length(ρ) > 100
    @test minimum(ρ) > 0.5
    @test maximum(ρ) < 6.0
    # hash documented
    readme = read(joinpath(@__DIR__, "data", "README.md"), String)
    @test occursin("SHA-256", readme)
end

@testset "absolute scoring" begin
    cases = Any[
        Dict(
            "case_type" => "smooth_order",
            "order_pass" => true,
            "diverged" => false,
            "nan_detected" => false,
            "positivity_ok" => true,
        ),
        Dict(
            "case_type" => "discontinuous",
            "excess_dissipation" => 0.05,
            "shock_thickness" => 4.0,
            "overshoot" => 0.05,
            "diverged" => false,
            "nan_detected" => false,
            "positivity_ok" => true,
        ),
    ]
    s = score_suite_absolute(cases)
    @test s["order_preservation"] == 1.0
    @test s["robustness"] == 1.0
    @test 0 < s["composite"] <= 1
    @test s["dissipation"] > 0.5
end

