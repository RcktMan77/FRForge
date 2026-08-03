using Test
using FRForge
using JSON

@testset "report schema skeleton" begin
    report = FRForge.report_skeleton()
    errs = validate_report_keys(report)
    @test isempty(errs)

    @test report["schema_version"] == 1
    @test report["package"] == "FRForge"
    @test report["method_name"] == "null"
    @test report["baseline_name"] === nothing
    @test report["overall_pass"] === true
    @test report["cases"] == []
    @test haskey(report["summary"], "scores")
    @test haskey(report["summary"]["scores"], "composite")
    @test !haskey(report, "overall_score")

    weights = report["scoring_weights"]
    @test weights["order_preservation"] ≈ 0.30
    @test weights["dissipation"] ≈ 0.25
    @test weights["shock_quality"] ≈ 0.25
    @test weights["robustness"] ≈ 0.20
    s = sum(values(weights))
    @test s ≈ 1.0
end

@testset "write and load report" begin
    mktempdir() do dir
        path = joinpath(dir, "nested", "report.json")
        report = write_report_skeleton(path; suite = "smoke")
        @test isfile(path)
        loaded = load_report(path)
        @test loaded["suite"] == "smoke"
        @test isempty(validate_report_keys(loaded))
        # Ensure JSON round-trip preserves schema version
        @test loaded["schema_version"] == 1
    end
end

@testset "validator catches missing keys" begin
    report = FRForge.report_skeleton()
    delete!(report, "schema_version")
    errs = validate_report_keys(report)
    @test any(contains(e, "schema_version") for e in errs)

    report2 = FRForge.report_skeleton()
    report2["overall_score"] = 0.5
    errs2 = validate_report_keys(report2)
    @test any(contains(e, "overall_score") for e in errs2)
end

@testset "cli test command" begin
    mktempdir() do dir
        path = joinpath(dir, "cli_report.json")
        # Simulate: frforge test --report path --suite smoke
        code = FRForge.main_cli(["test", "--report", path, "--suite", "smoke"])
        @test code == 0
        @test isfile(path)
        loaded = load_report(path)
        @test loaded["command"] == "test"
        @test loaded["suite"] == "smoke"
        @test isempty(validate_report_keys(loaded))
    end
end
