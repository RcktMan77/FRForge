using Test
using FRForge
using Dates

@testset "experiment log paths and frozen scheme" begin
    @test FROZEN_INVENT_SCHEME.points == "GL"
    @test FROZEN_INVENT_SCHEME.flux == "Rusanov"
    @test FROZEN_INVENT_SCHEME.time == "SSP-RK3"
    @test isfile(default_experiment_log_path())
    @test isdir(joinpath(package_root(), "research"))
    ids = list_experiment_entry_ids(default_experiment_log_path())
    @test "20260803-persson_av-baseline" in ids
    @test "20260803-scaled_persson-invent" in ids
    @test narrative_required("promising")
    @test narrative_required("accepted_candidate")
    @test !narrative_required("pass_gates")
    @test !narrative_required("rejected")
end

@testset "entry format and append" begin
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
    met = deepcopy(bas)
    met["method_name"] = "toy_method"
    met["baseline_name"] = "persson_av"
    met["cases"][2]["name"] = "sod_p2_ne64_toy"
    met["cases"][2]["excess_dissipation"] = 0.01
    met["cases"][2]["overshoot"] = 0.01
    met["cases"][2]["shock_thickness"] = 3.0
    met["summary"]["scores"] = score_suite_absolute(met["cases"])
    cmp = classify_candidate(met, bas; δ=0.02, vtk_produced=false)

    entry = entry_from_invent(
        "toy_method",
        met,
        bas,
        cmp;
        hypothesis="Test hypothesis for log unit test.",
        lessons="Test lessons for log unit test.",
        artifacts=Dict("method_report" => "tmp/method.json"),
        entry_id="testdate-toy_method-unit",
        date=Date(2026, 8, 3),
    )
    @test entry["id"] == "testdate-toy_method-unit"
    @test entry["scheme"]["points"] == "GL"
    @test entry["metrics"]["candidate_status"] == cmp["candidate_status"]
    md = format_entry_markdown(entry)
    @test occursin("### testdate-toy_method-unit", md)
    @test occursin("Test hypothesis", md)
    @test occursin("points=GL", md)

    mktempdir() do dir
        logp = joinpath(dir, "experiment_log.md")
        write(logp, "# Test log\n\n## Entries\n\n")
        ymlp = joinpath(dir, "experiment_log.yaml")
        write(ymlp, "entries:\n")
        append_experiment_entry!(logp, entry; yaml_path=ymlp)
        ids = list_experiment_entry_ids(logp)
        @test "testdate-toy_method-unit" in ids
        body = read(logp, String)
        @test occursin("Test lessons", body)
        ybody = read(ymlp, String)
        @test occursin("testdate-toy_method-unit", ybody)

        # promising without narrative → placeholders
        cmp_p = deepcopy(cmp)
        # force promising-like path if not already
        if cmp_p["candidate_status"] != "promising"
            # still test entry_from_invent with fake promising status
            cmp_p["candidate_status"] = "promising"
        end
        e2 = entry_from_invent("toy_method", met, bas, cmp_p; entry_id="x-promising")
        @test e2["narrative_complete"] == false || (
            e2["hypothesis"] != FRForge.NARRATIVE_PLACEHOLDER &&
            e2["lessons"] != FRForge.NARRATIVE_PLACEHOLDER
        )
        e3 = entry_from_invent(
            "toy_method",
            met,
            bas,
            cmp_p;
            entry_id="x-promising-full",
            hypothesis="",
            lessons="",
        )
        # status forced promising → placeholders
        cmp_p2 = Dict{String,Any}(cmp)
        cmp_p2["candidate_status"] = "promising"
        e4 = entry_from_invent("toy_method", met, bas, cmp_p2; entry_id="x4")
        @test e4["hypothesis"] == FRForge.NARRATIVE_PLACEHOLDER
        @test e4["lessons"] == FRForge.NARRATIVE_PLACEHOLDER
        @test e4["narrative_complete"] == false
    end
end

@testset "cli log list" begin
    code = main_cli(["log", "list"])
    @test code == 0
end
