# Fine-mesh confirm: presets, classify, log entry (no heavy suite in required CI).

@testset "confirm mesh presets" begin
    for name in ("confirm", "presentation", "quick")
        spec = get_confirm_preset(name)
        @test spec.riemann_n >= 8
        @test spec.dmr_nx >= 8
        note = mesh_note_string(spec, name)
        @test occursin("riemann", note)
        @test occursin("dmr", note)
        @test occursin("preset=$name", note)
    end
    @test_throws ErrorException get_confirm_preset("nope")
    sumd = mesh_summary_dict(get_confirm_preset("quick"), "quick"; include_smooth=false)
    @test !haskey(sumd, "vortex")
    @test sumd["riemann"]["n"] == 16
end

@testset "classify_confirm hard gates" begin
    function mk_report(;
        method="m",
        overall=true,
        diverged=false,
        nan=false,
        r_pass=true,
        d_pass=true,
        order_pass_flag=true,
        S_order=1.0,
        S_rob=1.0,
    )
        cases = Any[
            Dict{String,Any}(
                "name" => "vortex",
                "case_type" => "smooth_order",
                "pass" => order_pass_flag,
                "order_pass" => order_pass_flag,
                "diverged" => false,
                "nan_detected" => false,
                "metrics" => Dict{String,Any}("confirm_role" => "smooth_order_core"),
            ),
            Dict{String,Any}(
                "name" => "riemann",
                "case_type" => "discontinuous",
                "pass" => r_pass,
                "diverged" => false,
                "nan_detected" => false,
                "metrics" => Dict{String,Any}("confirm_role" => "riemann_cfg6"),
            ),
            Dict{String,Any}(
                "name" => "dmr",
                "case_type" => "discontinuous",
                "pass" => d_pass,
                "diverged" => false,
                "nan_detected" => false,
                "metrics" => Dict{String,Any}("confirm_role" => "double_mach_reduced"),
            ),
        ]
        return Dict{String,Any}(
            "method_name" => method,
            "overall_pass" => overall,
            "diverged" => diverged,
            "nan_detected" => nan,
            "hard_gate_failures" => String[],
            "cases" => cases,
            "summary" => Dict{String,Any}(
                "scores" => Dict{String,Any}(
                    "order_preservation" => S_order,
                    "dissipation" => 0.8,
                    "shock_quality" => 0.8,
                    "robustness" => S_rob,
                    "composite" => 0.85,
                ),
            ),
            "mesh_summary" => Dict{String,Any}("preset" => "quick"),
        )
    end

    bas = mk_report(; method="persson_av")
    met = mk_report(; method="scaled_persson")
    cmp = classify_confirm(met, bas; preset="quick")
    @test cmp["confirmation_status"] == "confirmed"
    @test all(r -> r["ok"] === true, cmp["rules"])

    # Baseline must finish — hard fail
    bas_bad = mk_report(; method="persson_av", overall=false, diverged=true)
    cmp2 = classify_confirm(met, bas_bad; preset="quick")
    @test cmp2["confirmation_status"] == "confirmation_failed"
    @test any(r -> r["id"] == "baseline_finished" && r["ok"] == false, cmp2["rules"])

    # Method multi-D fail
    met_bad = mk_report(; method="scaled_persson", r_pass=false, overall=false)
    cmp3 = classify_confirm(met_bad, bas; preset="quick")
    @test cmp3["confirmation_status"] == "confirmation_failed"

    # Order regression
    met_ord = mk_report(; method="scaled_persson", S_order=0.0, order_pass_flag=false)
    bas_ord = mk_report(; method="persson_av", S_order=1.0)
    # method overall still true but order rule fails
    met_ord["overall_pass"] = true
    cmp4 = classify_confirm(met_ord, bas_ord; preset="quick")
    @test cmp4["confirmation_status"] == "confirmation_failed"
    @test any(r -> r["id"] == "order_no_regression" && r["ok"] == false, cmp4["rules"])
end

@testset "entry_from_confirm log shape" begin
    bas = Dict{String,Any}(
        "method_name" => "persson_av",
        "overall_pass" => true,
        "diverged" => false,
        "nan_detected" => false,
        "cases" => Any[],
        "summary" => Dict("scores" => Dict("order_preservation" => 1.0, "robustness" => 1.0)),
    )
    met = deepcopy(bas)
    met["method_name"] = "scaled_persson"
    met["mesh_summary"] = mesh_summary_dict(get_confirm_preset("quick"), "quick")
    cmp = Dict{String,Any}(
        "confirmation_status" => "confirmed",
        "baseline_name" => "persson_av",
        "method_overall_pass" => true,
        "baseline_overall_pass" => true,
        "absolute_scores" => Dict("order_preservation" => 1.0, "robustness" => 1.0, "composite" => 0.9),
        "mesh_summary" => met["mesh_summary"],
        "rules" => Any[],
    )
    entry = entry_from_confirm("scaled_persson", met, bas, cmp; preset="quick")
    @test entry["status"] == "confirmed"
    @test entry["metrics"]["confirmation_status"] == "confirmed"
    @test occursin("riemann", entry["metrics"]["mesh"])
    md = format_entry_markdown(entry)
    @test occursin("confirmation_status", md)
    @test occursin("confirmed", md)
    @test occursin("mesh:", md)
end

@testset "log summary surfaces confirmation_failed" begin
    text = """
### 20260804-demo-confirm

- **date:** 2026-08-04
- **method:** demo_method
- **baseline:** persson_av
- **hypothesis:** test
- **scheme:** points=GL, flux=Rusanov, time=SSP-RK3
- **metrics:**
  - confirmation_status: confirmation_failed
  - mesh: preset=quick; riemann 16²
- **lessons:** blew up
- **status:** confirmation_failed
- **artifacts:** none

"""
    entries = parse_experiment_log_text(text)
    s = log_summary(entries)
    @test haskey(s, "confirmation_failed")
    @test length(s["confirmation_failed"]) >= 1
    out = format_log_summary_text(s)
    @test occursin("confirmation_failed", out)
    @test occursin("CAUTION", out)
end

@testset "method_has_confirm_pass" begin
    mktempdir() do dir
        logp = joinpath(dir, "log.md")
        write(logp, """
### 20260804-m-confirm

- **date:** 2026-08-04
- **method:** my_method
- **baseline:** persson_av
- **hypothesis:** h
- **scheme:** points=GL, flux=Rusanov, time=SSP-RK3
- **metrics:**
  - confirmation_status: confirmed
- **lessons:** ok
- **status:** confirmed
- **artifacts:** none
""")
        ok, note = method_has_confirm_pass("my_method"; log_path=logp)
        @test ok
        @test occursin("confirmed", note)
        ok2, _ = method_has_confirm_pass("other"; log_path=logp)
        @test !ok2
    end
end
