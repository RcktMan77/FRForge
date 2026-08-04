# Robustness matrix: re-evaluate capturing methods across scheme axes.
#
# CI policy: `matrix=:ci` uses a reduced light suite + small cell set.
# Full product is local / nightly only (`matrix=:full`).

const ROBUSTNESS_POINTS = (:gl, :gll)
const ROBUSTNESS_FLUXES = (:rusanov, :hllc)
const ROBUSTNESS_TIMES = (:ssp_rk3, :ssp_rk2)

"""Statuses considered publication-grade on the default invent scheme."""
const PROMOTION_INVENT_STATUSES = Set(["promising", "accepted_candidate"])

"""
    scheme_slug(scheme) -> String

Filesystem-safe slug: `gl_rusanov_ssp_rk3`.
"""
function scheme_slug(scheme::SchemeConfig)
    return string(scheme.points, "_", scheme.flux, "_", scheme.time)
end

"""
    robustness_cells(matrix=:full) -> Vector{SchemeConfig}

- `:full` — product of GL/GLL × Rusanov/HLLC × SSP-RK3/SSP-RK2 (8 cells)
- `:ci` — default scheme + one hard corner (GLL × HLLC × SSP-RK3)
"""
function robustness_cells(matrix::Symbol=:full)
    if matrix === :full
        cells = SchemeConfig[]
        for p in ROBUSTNESS_POINTS, f in ROBUSTNESS_FLUXES, t in ROBUSTNESS_TIMES
            push!(cells, SchemeConfig(; points=p, flux=f, time=t))
        end
        return cells
    elseif matrix === :ci
        return [
            DEFAULT_SCHEME,
            SchemeConfig(; points=:gll, flux=:hllc, time=:ssp_rk3),
        ]
    else
        error("Unknown matrix kind $matrix (use :full or :ci)")
    end
end

"""
    cell_ok(report) -> Bool

True if the cell did not diverge / NaN and overall_pass.
"""
function cell_ok(report::AbstractDict)
    return get(report, "overall_pass", false) === true &&
           get(report, "diverged", false) !== true &&
           get(report, "nan_detected", false) !== true
end

"""
    order_preserved(report) -> Bool

True if every smooth_order case in the report has order_pass.
"""
function order_preserved(report::AbstractDict)
    cases = get(report, "cases", Any[])
    smooth = filter(c -> get(c, "case_type", "") == "smooth_order", cases)
    isempty(smooth) && return true
    return all(c -> get(c, "order_pass", false) === true, smooth)
end

"""
    evaluate_robustness_cell(method_name, baseline_name, scheme; light, δ) -> Dict

Run quant (or light) suite for method and baseline under `scheme`, classify, return cell dict.
"""
function evaluate_robustness_cell(
    method_name::AbstractString,
    baseline_name::AbstractString,
    scheme::SchemeConfig;
    light::Bool=false,
    δ::Real=DEFAULT_SCORE_MARGIN,
    vtk_produced::Bool=false,
)
    suite = light ? :light : :quant
    t0 = time()
    met = run_method_report(method_name; suite=suite, scheme=scheme)
    bas = run_method_report(baseline_name; suite=suite, scheme=scheme)
    cmp = classify_candidate(met, bas; δ=δ, vtk_produced=vtk_produced)
    ok = cell_ok(met)
    ord = order_preserved(met)
    cell_status = if !ok
        "failed"
    elseif !ord
        "order_fail"
    else
        "ok"
    end
    scores = get(get(met, "summary", Dict()), "scores", Dict())
    return Dict{String,Any}(
        "scheme" => scheme_dict(scheme),
        "scheme_slug" => scheme_slug(scheme),
        "cell_status" => cell_status,
        "ok" => ok && ord,
        "overall_pass" => get(met, "overall_pass", false),
        "diverged" => get(met, "diverged", false),
        "nan_detected" => get(met, "nan_detected", false),
        "order_preserved" => ord,
        "candidate_status" => cmp["candidate_status"],
        "composite" => get(cmp, "composite", nothing),
        "baseline_composite" => get(cmp, "baseline_composite", nothing),
        "composite_margin" => get(cmp, "composite_margin", nothing),
        "tradeoff_ok" => get(cmp, "tradeoff_ok", nothing),
        "absolute_scores" => scores,
        "wall_time_sec" => time() - t0,
        "method_report" => met,
        "baseline_report" => bas,
        "comparison" => cmp,
    )
end

"""
    assess_publication_grade(cells; narrative_complete=false, invent_status=nothing,
                             confirm_passed=false) -> Dict

Promotion rule (all must hold for `eligible=true`):

1. **Default scheme** cell is `ok` and invent-status is `promising` or better
   (either from that cell's `candidate_status` or explicit `invent_status`).
2. **HLLC cells** (less-dissipative flux): all `ok` with order preserved.
3. **GLL cells**: all `ok` (no catastrophic failure).
4. **`narrative_complete`**: hypothesis + lessons filled (caller responsibility).
5. **`confirm_passed`**: fine-mesh `frforge confirm` status is `confirmed`
   (short-list on coarse invent alone is not enough for publication_grade).

Does **not** auto-promote; returns an assessment for the log / agents.
"""
function assess_publication_grade(
    cells::AbstractVector;
    narrative_complete::Bool=false,
    invent_status::Union{Nothing,AbstractString}=nothing,
    confirm_passed::Bool=false,
)
    rules = Dict{String,Any}[]
    default_cells = filter(c -> _is_default_scheme(c), cells)
    hllc_cells = filter(c -> _scheme_flux(c) === :hllc, cells)
    gll_cells = filter(c -> _scheme_points(c) === :gll, cells)

    # Rule 1: default scheme
    r1_ok = false
    r1_notes = "no default-scheme cell in matrix"
    if !isempty(default_cells)
        dc = default_cells[1]
        status = something(invent_status, string(get(dc, "candidate_status", "")))
        r1_ok =
            get(dc, "ok", false) === true &&
            (status in PROMOTION_INVENT_STATUSES || status == "publication_grade")
        r1_notes =
            "default cell_status=$(get(dc, "cell_status", "?")) invent_status=$status ok=$(get(dc, "ok", false))"
    end
    push!(
        rules,
        Dict(
            "id" => "default_promising",
            "ok" => r1_ok,
            "notes" => r1_notes,
        ),
    )

    # Rule 2: HLLC
    r2_ok = isempty(hllc_cells) ? false : all(c -> get(c, "ok", false) === true, hllc_cells)
    r2_notes = isempty(hllc_cells) ? "no HLLC cells" :
               "hllc_ok=$(count(c -> get(c, "ok", false) === true, hllc_cells))/$(length(hllc_cells))"
    push!(rules, Dict("id" => "hllc_robust", "ok" => r2_ok, "notes" => r2_notes))

    # Rule 3: GLL
    r3_ok = isempty(gll_cells) ? false : all(c -> get(c, "ok", false) === true, gll_cells)
    r3_notes = isempty(gll_cells) ? "no GLL cells" :
               "gll_ok=$(count(c -> get(c, "ok", false) === true, gll_cells))/$(length(gll_cells))"
    push!(rules, Dict("id" => "gll_stable", "ok" => r3_ok, "notes" => r3_notes))

    # Rule 4: narrative
    push!(
        rules,
        Dict(
            "id" => "narrative_complete",
            "ok" => narrative_complete,
            "notes" => narrative_complete ? "hypothesis+lessons provided" : "hypothesis/lessons incomplete",
        ),
    )

    # Rule 5: fine-mesh confirm
    push!(
        rules,
        Dict(
            "id" => "fine_mesh_confirm",
            "ok" => confirm_passed,
            "notes" => confirm_passed ?
                       "frforge confirm passed (confirmed)" :
                       "missing fine-mesh confirm — run frforge confirm after short-list",
        ),
    )

    all_ok = all(r -> r["ok"] === true, rules)
    return Dict{String,Any}(
        "eligible" => all_ok,
        "recommended_status" => all_ok ? "publication_grade" : "robustness_pending",
        "rules" => rules,
        "n_cells" => length(cells),
        "n_ok" => count(c -> get(c, "ok", false) === true, cells),
    )
end

function _scheme_field_symbol(cell::AbstractDict, human_key::String, sym_key::String)
    sch = get(cell, "scheme", Dict())
    if haskey(sch, sym_key)
        return Symbol(lowercase(string(sch[sym_key])))
    end
    raw = lowercase(replace(string(get(sch, human_key, "")), "-" => "_"))
    # map human labels
    raw == "gl" && return :gl
    raw == "gll" && return :gll
    raw == "rusanov" && return :rusanov
    raw == "hllc" && return :hllc
    raw == "ssp_rk3" && return :ssp_rk3
    raw == "ssp_rk2" && return :ssp_rk2
    return Symbol(raw)
end

_scheme_points(cell::AbstractDict) = _scheme_field_symbol(cell, "points", "points_symbol")
_scheme_flux(cell::AbstractDict) = _scheme_field_symbol(cell, "flux", "flux_symbol")
_scheme_time(cell::AbstractDict) = _scheme_field_symbol(cell, "time", "time_symbol")

function _is_default_scheme(cell::AbstractDict)
    return _scheme_points(cell) === :gl &&
           _scheme_flux(cell) === :rusanov &&
           _scheme_time(cell) === :ssp_rk3
end

"""
    run_robustness_matrix(method_name; baseline, matrix, report_dir, light, ...) -> Dict

Evaluate all cells, write per-cell JSON + summary, optionally append experiment log.
"""
function run_robustness_matrix(
    method_name::AbstractString;
    baseline::AbstractString="persson_av",
    matrix::Symbol=:ci,
    report_dir::AbstractString="results/robustness",
    light::Union{Nothing,Bool}=nothing,
    δ::Real=DEFAULT_SCORE_MARGIN,
    append_log::Bool=true,
    log_path::Union{Nothing,AbstractString}=nothing,
    narrative_complete::Bool=false,
    invent_status::Union{Nothing,AbstractString}=nothing,
    hypothesis::AbstractString="",
    lessons::AbstractString="",
)
    # CI matrix always light; full matrix defaults to full quant unless light=true
    use_light = light === nothing ? (matrix === :ci) : light
    cells_cfg = robustness_cells(matrix)
    out_root = joinpath(report_dir, String(method_name))
    mkpath(out_root)

    println("Robustness matrix: method=$method_name baseline=$baseline matrix=$matrix light=$use_light")
    println("  cells: $(length(cells_cfg))  →  $out_root")

    cell_results = Dict{String,Any}[]
    for (i, scheme) in enumerate(cells_cfg)
        slug = scheme_slug(scheme)
        println("  [$i/$(length(cells_cfg))] $slug ...")
        cell = evaluate_robustness_cell(
            method_name,
            baseline,
            scheme;
            light=use_light,
            δ=δ,
        )
        # Strip bulky nested reports from summary cells; write full reports separately
        met = cell["method_report"]
        bas = cell["baseline_report"]
        cmp = cell["comparison"]
        # Per-cell slug is unique; keep historical filenames (not invent trio layout).
        met_path = joinpath(out_root, "method_$(slug).json")
        bas_path = joinpath(out_root, "baseline_$(slug).json")
        cmp_path = joinpath(out_root, "compare_$(slug).json")
        write_report(met_path, met)
        write_report(bas_path, bas)
        write_json_pretty(cmp_path, cmp)
        slim = Dict{String,Any}(
            "scheme" => cell["scheme"],
            "scheme_slug" => slug,
            "cell_status" => cell["cell_status"],
            "ok" => cell["ok"],
            "overall_pass" => cell["overall_pass"],
            "diverged" => cell["diverged"],
            "nan_detected" => cell["nan_detected"],
            "order_preserved" => cell["order_preserved"],
            "candidate_status" => cell["candidate_status"],
            "composite" => cell["composite"],
            "baseline_composite" => cell["baseline_composite"],
            "composite_margin" => cell["composite_margin"],
            "tradeoff_ok" => cell["tradeoff_ok"],
            "absolute_scores" => cell["absolute_scores"],
            "wall_time_sec" => cell["wall_time_sec"],
            "artifacts" => Dict(
                "method_report" => met_path,
                "baseline_report" => bas_path,
                "compare" => cmp_path,
            ),
        )
        push!(cell_results, slim)
        println(
            "    status=$(slim["cell_status"]) composite=$(slim["composite"]) invent=$(slim["candidate_status"])",
        )
    end

    conf_ok, conf_note = method_has_confirm_pass(method_name)
    promotion = assess_publication_grade(
        cell_results;
        narrative_complete=narrative_complete,
        invent_status=invent_status,
        confirm_passed=conf_ok,
    )
    # attach note for agents
    for r in promotion["rules"]
        if r["id"] == "fine_mesh_confirm"
            r["notes"] = conf_note
            r["ok"] = conf_ok
        end
    end
    promotion["eligible"] = all(r -> r["ok"] === true, promotion["rules"])
    promotion["recommended_status"] =
        promotion["eligible"] ? "publication_grade" : "robustness_pending"

    summary = Dict{String,Any}(
        "schema_version" => SCHEMA_VERSION,
        "command" => "robustness",
        "method_name" => String(method_name),
        "baseline_name" => String(baseline),
        "matrix" => string(matrix),
        "light" => use_light,
        "timestamp_utc" =>
            Dates.format(Dates.now(Dates.UTC), dateformat"yyyy-mm-ddTHH:MM:SS") * "Z",
        "git_commit" => git_commit_short(),
        "cells" => cell_results,
        "promotion" => promotion,
        "score_margin_threshold" => float(δ),
        "ci_policy" => "required CI uses matrix=ci + light suite only; full matrix is local/nightly",
    )

    summary_path = joinpath(out_root, "summary.json")
    write_json_pretty(summary_path, summary)
    println("Summary: $summary_path")
    println(
        "Promotion: eligible=$(promotion["eligible"]) → $(promotion["recommended_status"])",
    )
    for r in promotion["rules"]
        println("  [$(r["ok"] ? "OK" : "FAIL")] $(r["id"]): $(r["notes"])")
    end

    if append_log
        lp = something(log_path, default_experiment_log_path())
        append_robustness_log_entry!(
            summary;
            log_path=lp,
            hypothesis=hypothesis,
            lessons=lessons,
        )
        println("Experiment log updated: $lp")
    end

    return summary
end

"""Append a Markdown robustness batch entry to the experiment log."""
function append_robustness_log_entry!(
    summary::AbstractDict;
    log_path::AbstractString=default_experiment_log_path(),
    hypothesis::AbstractString="",
    lessons::AbstractString="",
)
    method = string(summary["method_name"])
    date = Dates.today()
    id = make_entry_id(method; date=date, suffix="robustness_$(summary["matrix"])")
    prom = summary["promotion"]
    hyp = isempty(strip(hypothesis)) ?
          "Robustness matrix ($(summary["matrix"]), light=$(summary["light"])) for method $method." :
          String(hypothesis)
    les = isempty(strip(lessons)) ?
          "See promotion rules in summary; eligible=$(prom["eligible"])." :
          String(lessons)

    lines = String[]
    push!(lines, "### $id")
    push!(lines, "")
    push!(lines, "- **date:** $date")
    push!(lines, "- **method:** $method")
    push!(lines, "- **baseline:** $(summary["baseline_name"])")
    push!(lines, "- **hypothesis:** $hyp")
    push!(lines, "- **scheme:** robustness matrix ($(summary["matrix"]))")
    push!(lines, "- **metrics:**")
    push!(lines, "  - matrix: $(summary["matrix"])")
    push!(lines, "  - light: $(summary["light"])")
    push!(lines, "  - n_cells: $(prom["n_cells"])")
    push!(lines, "  - n_ok: $(prom["n_ok"])")
    push!(lines, "  - promotion_eligible: $(prom["eligible"])")
    push!(lines, "  - recommended_status: $(prom["recommended_status"])")
    for c in summary["cells"]
        push!(
            lines,
            "  - cell $(c["scheme_slug"]): $(c["cell_status"]) composite=$(c["composite"]) invent=$(c["candidate_status"])",
        )
    end
    for r in prom["rules"]
        push!(lines, "  - rule $(r["id"]): $(r["ok"]) ($(r["notes"]))")
    end
    push!(lines, "- **lessons:** $les")
    push!(
        lines,
        "- **status:** $(prom["eligible"] ? "publication_grade" : "robustness_pending")",
    )
    push!(lines, "- **artifacts:**")
    push!(
        lines,
        "  - summary: results/robustness/$(method)/summary.json",
    )
    push!(lines, "")

    mkpath(dirname(log_path))
    open(log_path, "a") do io
        write(io, "\n")
        write(io, join(lines, "\n"))
        write(io, "\n")
    end
    return id
end
