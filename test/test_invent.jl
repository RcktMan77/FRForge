using Test
using FRForge

@testset "method registry lists invent methods" begin
    names = list_methods()
    @test "null" in names
    @test "persson_av" in names
    @test "scaled_persson" in names
    m = get_capturing_method("scaled_persson")
    @test m isa ScaledPerssonMethod
    @test method_params(m)["c_av"] > method_params(PerssonAVMethod())["c_av"]
end

@testset "classify_candidate statuses" begin
    bas = FRForge.report_skeleton(;
        command="invent",
        suite="quant",
        method_name="persson_av",
        overall_pass=true,
        cases=Any[
            Dict(
                "name" => "euler_smooth_order_p2",
                "case_type" => "smooth_order",
                "order_pass" => true,
                "diverged" => false,
                "nan_detected" => false,
                "positivity_ok" => true,
                "equation" => "euler1d",
                "p" => 2,
                "capturing_method" => "null",
                "pass" => true,
                "conservation_residual" => 0.0,
                "conservation_pass" => true,
                "conservation_metric" => "periodic_mass_change",
                "wall_time_sec" => 0.0,
                "metrics" => Dict{String,Any}(),
            ),
            Dict(
                "name" => "sod_p2_ne64_persson_av",
                "case_type" => "discontinuous",
                "excess_dissipation" => 0.05,
                "shock_thickness" => 5.0,
                "overshoot" => 0.05,
                "diverged" => false,
                "nan_detected" => false,
                "positivity_ok" => true,
                "equation" => "euler1d",
                "p" => 2,
                "capturing_method" => "persson_av",
                "pass" => true,
                "conservation_residual" => 0.0,
                "conservation_pass" => true,
                "conservation_metric" => "none",
                "wall_time_sec" => 0.0,
                "metrics" => Dict{String,Any}(),
            ),
        ],
        fill_scores=true,
    )
    # Better method: lower excess dissip and overshoot
    met = deepcopy(bas)
    met["method_name"] = "scaled_persson"
    met["cases"][2]["name"] = "sod_p2_ne64_scaled_persson"
    met["cases"][2]["excess_dissipation"] = 0.01
    met["cases"][2]["overshoot"] = 0.01
    met["cases"][2]["shock_thickness"] = 3.0
    met["summary"]["scores"] = score_suite_absolute(met["cases"])

    cmp = classify_candidate(met, bas; δ=0.02, vtk_produced=false)
    @test cmp["candidate_status"] in ("promising", "pass_gates", "accepted_candidate")
    @test cmp["composite_margin"] > 0
    @test cmp["tradeoff_ok"]

    cmp2 = classify_candidate(met, bas; δ=0.02, vtk_produced=true)
    if cmp["candidate_status"] == "promising"
        @test cmp2["candidate_status"] == "accepted_candidate"
    end

    # Rejected if gates fail
    bad = deepcopy(met)
    bad["overall_pass"] = false
    bad["diverged"] = true
    cmp3 = classify_candidate(bad, bas)
    @test cmp3["candidate_status"] == "rejected"
end

@testset "score_reports round-trip" begin
    mktempdir() do dir
        bas = FRForge.report_skeleton(;
            method_name="persson_av",
            overall_pass=true,
            cases=Any[
                Dict(
                    "name" => "sod_base",
                    "case_type" => "discontinuous",
                    "excess_dissipation" => 0.1,
                    "shock_thickness" => 6.0,
                    "overshoot" => 0.08,
                    "diverged" => false,
                    "nan_detected" => false,
                    "positivity_ok" => true,
                    "equation" => "euler1d",
                    "p" => 2,
                    "capturing_method" => "persson_av",
                    "pass" => true,
                    "conservation_residual" => 0.0,
                    "conservation_pass" => true,
                    "conservation_metric" => "none",
                    "wall_time_sec" => 0.0,
                    "metrics" => Dict{String,Any}(),
                ),
            ],
            fill_scores=true,
        )
        met = deepcopy(bas)
        met["method_name"] = "scaled_persson"
        met["cases"][1]["excess_dissipation"] = 0.02
        met["cases"][1]["overshoot"] = 0.02
        met["summary"]["scores"] = score_suite_absolute(met["cases"])
        bp = joinpath(dir, "b.json")
        mp = joinpath(dir, "m.json")
        write_report(bp, bas)
        write_report(mp, met)
        cmp = score_reports(mp, bp; δ=0.01)
        @test haskey(cmp, "candidate_status")
        @test cmp["composite"] > cmp["baseline_composite"]
    end
end

