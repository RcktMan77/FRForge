using Test
using FRForge
using Dates

@testset "robustness cell enumeration" begin
    full = robustness_cells(:full)
    @test length(full) == 8
    @test DEFAULT_SCHEME in full
    @test any(s -> s.points === :gll && s.flux === :hllc && s.time === :ssp_rk3, full)
    ci = robustness_cells(:ci)
    @test length(ci) == 2
    @test ci[1] == DEFAULT_SCHEME
    @test scheme_slug(DEFAULT_SCHEME) == "gl_rusanov_ssp_rk3"
end

@testset "assess_publication_grade rules" begin
    function fake_cell(; points, flux, time, ok, cand = "pass_gates")
        sch = SchemeConfig(; points = points, flux = flux, time = time)
        return Dict{String,Any}(
            "scheme" => scheme_dict(sch),
            "scheme_slug" => scheme_slug(sch),
            "ok" => ok,
            "cell_status" => ok ? "ok" : "failed",
            "candidate_status" => cand,
        )
    end
    # All ok + promising default + narrative
    cells = [
        fake_cell(; points = :gl, flux = :rusanov, time = :ssp_rk3, ok = true, cand = "promising"),
        fake_cell(; points = :gll, flux = :hllc, time = :ssp_rk3, ok = true, cand = "pass_gates"),
        fake_cell(; points = :gl, flux = :hllc, time = :ssp_rk2, ok = true, cand = "pass_gates"),
        fake_cell(;
            points = :gll,
            flux = :rusanov,
            time = :ssp_rk3,
            ok = true,
            cand = "pass_gates",
        ),
    ]
    a = assess_publication_grade(cells; narrative_complete = true, confirm_passed = true)
    @test a["eligible"] == true
    @test a["recommended_status"] == "publication_grade"

    # Missing fine-mesh confirm
    a0 = assess_publication_grade(cells; narrative_complete = true, confirm_passed = false)
    @test a0["eligible"] == false
    @test any(r -> r["id"] == "fine_mesh_confirm" && r["ok"] == false, a0["rules"])

    # HLLC failure
    cells_bad = copy(cells)
    cells_bad[2] = fake_cell(; points = :gll, flux = :hllc, time = :ssp_rk3, ok = false)
    a2 = assess_publication_grade(cells_bad; narrative_complete = true, confirm_passed = true)
    @test a2["eligible"] == false
    @test any(r -> r["id"] == "hllc_robust" && r["ok"] == false, a2["rules"])

    # Default not promising
    cells3 = [
        fake_cell(; points = :gl, flux = :rusanov, time = :ssp_rk3, ok = true, cand = "pass_gates"),
        fake_cell(; points = :gll, flux = :hllc, time = :ssp_rk3, ok = true),
    ]
    a3 = assess_publication_grade(cells3; narrative_complete = true, confirm_passed = true)
    @test a3["eligible"] == false

    # Narrative incomplete
    a4 = assess_publication_grade(cells; narrative_complete = false, confirm_passed = true)
    @test a4["eligible"] == false
end

@testset "order_preserved and cell_ok helpers" begin
    good = Dict(
        "overall_pass" => true,
        "diverged" => false,
        "nan_detected" => false,
        "cases" => Any[
            Dict("case_type" => "smooth_order", "order_pass" => true),
            Dict("case_type" => "discontinuous", "order_pass" => false),
        ],
    )
    @test cell_ok(good)
    @test order_preserved(good)
    bad = deepcopy(good)
    bad["diverged"] = true
    @test !cell_ok(bad)
    bad2 = deepcopy(good)
    bad2["cases"][1]["order_pass"] = false
    @test !order_preserved(bad2)
end

@testset "robustness light cell smoke (CI-tier)" begin
    # Single light cell under DEFAULT_SCHEME — keeps required CI budget.
    # Compare persson_av vs itself so classification is well-defined.
    cell = evaluate_robustness_cell("persson_av", "persson_av", DEFAULT_SCHEME; light = true)
    @test haskey(cell, "cell_status")
    @test haskey(cell, "scheme_slug")
    @test cell["scheme_slug"] == "gl_rusanov_ssp_rk3"
    @test cell["cell_status"] in ("ok", "order_fail", "failed")
    # Self-baseline: composite margin ~0
    @test cell["composite"] isa Real
end

@testset "run_robustness_matrix ci smoke" begin
    mktempdir() do dir
        logp = joinpath(dir, "experiment_log.md")
        write(logp, "# test log\n")
        # One-method matrix=ci is two light cells — acceptable for unit tests
        summary = run_robustness_matrix(
            "null";
            baseline = "null",
            matrix = :ci,
            report_dir = joinpath(dir, "rob"),
            light = true,
            append_log = true,
            log_path = logp,
            narrative_complete = false,
            invent_status = "pass_gates",
        )
        @test summary["matrix"] == "ci"
        @test length(summary["cells"]) == 2
        @test isfile(joinpath(dir, "rob", "null", "summary.json"))
        @test haskey(summary, "promotion")
        body = read(logp, String)
        @test occursin("robustness_ci", body)
    end
end
