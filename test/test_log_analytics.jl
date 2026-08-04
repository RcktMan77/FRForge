using Test
using FRForge
using Dates

@testset "parse seed experiment log" begin
    path = default_experiment_log_path()
    @test isfile(path)
    entries = parse_experiment_log(path)
    @test length(entries) >= 3
    ids = [e["id"] for e in entries]
    @test "20260803-persson_av-baseline" in ids
    @test "20260803-scaled_persson-invent" in ids
    @test "20260803-m0-m8-platform" in ids

    bas = get_experiment_entry(entries, "20260803-persson_av-baseline")
    @test bas !== nothing
    @test bas["method"] == "persson_av"
    @test bas["status"] == "baseline"
    @test get(bas["metrics"], "candidate_status", "") == "baseline"
    @test bas["scheme"]["points"] == "GL"
    @test bas["scheme"]["flux"] == "Rusanov"
    @test haskey(bas["artifacts"], "method_report")

    sc = get_experiment_entry(entries, "20260803-scaled_persson-invent")
    @test sc !== nothing
    @test sc["method"] == "scaled_persson"
    @test get(sc["metrics"], "candidate_status", "") == "pass_gates"
    @test sc["metrics"]["composite"] isa Number
    @test sc["metrics"]["composite"] > 0.9
    @test occursin("NullCapturing", sc["lessons"]) || occursin("structural", lowercase(sc["lessons"]))

    plat = get_experiment_entry(entries, "20260803-m0-m8-platform")
    @test plat !== nothing
    @test plat["kind"] == "platform"
end

@testset "parse messy log fragment" begin
    messy = """
# Preamble

### 20260101-platform-note

- **date:** 2026-01-01

- **method:** _(platform milestone, not a capturing method)_

- **baseline:** n/a
- **hypothesis:** Platform only.
- **scheme:** points = GL , flux=Rusanov, time = SSP-RK3
- **metrics:** n/a
- **status:** archived
- **artifacts:** none

### 20260102-toy_method-invent

- **date:** 2026-01-02
- *method*: toy_method
- **baseline:** persson_av
- **hypothesis:** Try a toy.
- **scheme:** points=GL, flux=Rusanov, time=SSP-RK3
- **metrics:**
  - candidate_status: pass_gates
  - composite: ~0.91 (approx)
  - order_preservation: 1.0
  - dissipation: 0.9
  - shock_quality: 0.7

- **lessons:** Keep c_av conservative.
- **status:** open
# no strengths / weaknesses on purpose
"""
    entries = parse_experiment_log_text(messy)
    @test length(entries) == 2
    plat = get_experiment_entry(entries, "20260101-platform-note")
    @test plat !== nothing
    @test plat["kind"] == "platform"
    @test plat["scheme"]["points"] == "GL"

    toy = get_experiment_entry(entries, "20260102-toy_method-invent")
    @test toy !== nothing
    @test toy["method"] == "toy_method"
    @test toy["metrics"]["composite"] isa Number
    @test isapprox(toy["metrics"]["composite"], 0.91; atol=1e-9)
    @test toy["metrics"]["candidate_status"] == "pass_gates"
    @test isempty(get(toy, "strengths", ""))  # missing optional
    @test occursin("c_av", toy["lessons"])
end

@testset "round-trip invent entry markdown" begin
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
    met["method_name"] = "rt_method"
    met["cases"][2]["excess_dissipation"] = 0.01
    met["summary"]["scores"] = score_suite_absolute(met["cases"])
    cmp = classify_candidate(met, bas; δ=0.02, vtk_produced=false)
    entry = entry_from_invent("rt_method", met, bas, cmp; hypothesis="H", lessons="L")
    md = format_entry_markdown(entry)
    parsed = parse_experiment_log_text("# Entries\n\n" * md)
    @test length(parsed) == 1
    p = parsed[1]
    @test p["method"] == "rt_method"
    @test p["metrics"]["composite"] isa Number
    @test isapprox(p["metrics"]["composite"], entry["metrics"]["composite"]; rtol=1e-12)
    @test p["hypothesis"] == "H"
    @test p["lessons"] == "L"
end

@testset "summary frontier pareto lessons views" begin
    entries = parse_experiment_log()
    summary = log_summary(entries)
    @test summary["n_entries"] >= 3
    @test haskey(summary["by_status"], "baseline")
    @test haskey(summary["latest_by_method"], "persson_av")
    @test haskey(summary["latest_by_method"], "scaled_persson")

    frontier = log_frontier(entries)
    methods = [r["method"] for r in frontier]
    @test "persson_av" in methods
    # scaled_persson is pass_gates with competitive composite → near-miss
    sc_rows = filter(r -> r["method"] == "scaled_persson", frontier)
    @test !isempty(sc_rows)
    @test sc_rows[1]["near"] === true || sc_rows[1]["candidate_status"] == "pass_gates"

    pareto = log_pareto(entries)
    # scaled_persson has full numeric scores
    @test any(r -> r["method"] == "scaled_persson", pareto)

    lessons = log_lessons(entries; query="NullCapturing")
    @test !isempty(lessons)
    @test any(r -> occursin("NullCapturing", r["text"]), lessons)

    lessons_all = log_lessons(entries)
    @test length(lessons_all) >= length(lessons)

    # formatters non-empty
    @test !isempty(format_log_summary_text(summary))
    @test !isempty(format_log_frontier_text(frontier))
    @test !isempty(format_log_lessons_text(lessons))

    e = get_experiment_entry(entries, "20260803-scaled_persson-invent")
    @test !isempty(format_log_entry_text(e))
end
