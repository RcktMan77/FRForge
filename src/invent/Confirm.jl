# Fine-mesh short-list confirmation (peer of invent / robustness).
#
# Coarse invent (Invent.jl) stays fast and owns composite-score history.
# `frforge confirm` re-runs key multi-D cases on finer meshes before
# publication-grade claims or snapshot freeze (--require-confirm).

"""Mesh / time parameters for one confirm case family."""
struct ConfirmMeshSpec
    riemann_n::Int
    riemann_p::Int
    riemann_t::Float64
    riemann_cfl::Float64
    dmr_nx::Int
    dmr_ny::Int
    dmr_p::Int
    dmr_t::Float64
    dmr_cfl::Float64
    dmr_Lx::Float64
    dmr_Ly::Float64
    vortex_p::Int
    vortex_n_list::Vector{Int}
    vortex_t::Float64
    vortex_cfl::Float64
    # Rough wall-time class on Apple Silicon single process (document; re-probe on change)
    wall_time_note::String
end

"""
Named mesh presets for confirmation.

- `confirm` (default): noticeably finer than CI-light; target ~10–30 min class.
- `presentation`: docs / paper-grade meshes (may take hours).
- `quick`: smoke / unit-test path.
"""
const CONFIRM_MESH_PRESETS = Dict{String,ConfirmMeshSpec}(
    "confirm" => ConfirmMeshSpec(
        64, 2, 0.15, 0.05,
        120, 40, 1, 0.08, 0.035,
        1.5, 0.5,
        2, [16, 32], 0.5, 0.08,
        "confirm preset: ~10–30 min typical (64² Riemann p=2 + DMR 120×40 + vortex); probe on your machine",
    ),
    "presentation" => ConfirmMeshSpec(
        192, 2, 0.15, 0.035,
        280, 100, 1, 0.08, 0.025,
        1.5, 0.5,
        2, [16, 32, 48], 0.5, 0.06,
        "presentation preset: paper-grade; may take hours (192² Riemann + DMR 280×100)",
    ),
    "quick" => ConfirmMeshSpec(
        16, 1, 0.06, 0.08,
        24, 8, 1, 0.03, 0.05,
        1.5, 0.5,
        2, [8, 16], 0.25, 0.1,
        "quick preset: smoke only (CI-scale meshes); not for publication claims",
    ),
)

function get_confirm_preset(name::AbstractString)
    key = lowercase(String(name))
    haskey(CONFIRM_MESH_PRESETS, key) ||
        error("Unknown confirm preset \"$name\" (use confirm|presentation|quick)")
    return CONFIRM_MESH_PRESETS[key]
end

function mesh_summary_dict(spec::ConfirmMeshSpec, preset::AbstractString; include_smooth::Bool=true)
    d = Dict{String,Any}(
        "preset" => String(preset),
        "riemann" => Dict(
            "n" => spec.riemann_n,
            "p" => spec.riemann_p,
            "t_final" => spec.riemann_t,
            "cfl" => spec.riemann_cfl,
        ),
        "dmr" => Dict(
            "nx" => spec.dmr_nx,
            "ny" => spec.dmr_ny,
            "p" => spec.dmr_p,
            "t_final" => spec.dmr_t,
            "cfl" => spec.dmr_cfl,
            "Lx" => spec.dmr_Lx,
            "Ly" => spec.dmr_Ly,
        ),
        "wall_time_note" => spec.wall_time_note,
    )
    if include_smooth
        d["vortex"] = Dict(
            "p" => spec.vortex_p,
            "n_list" => copy(spec.vortex_n_list),
            "t_final" => spec.vortex_t,
            "cfl" => spec.vortex_cfl,
        )
    end
    return d
end

function mesh_note_string(spec::ConfirmMeshSpec, preset::AbstractString; include_smooth::Bool=true)
    parts = [
        "preset=$preset",
        "riemann $(spec.riemann_n)² p=$(spec.riemann_p) t=$(spec.riemann_t)",
        "dmr $(spec.dmr_nx)×$(spec.dmr_ny) p=$(spec.dmr_p) t=$(spec.dmr_t)",
    ]
    if include_smooth
        push!(
            parts,
            "vortex p=$(spec.vortex_p) n=$(join(spec.vortex_n_list, ","))",
        )
    end
    return join(parts, "; ")
end

"""
    run_confirm_suite(method_name; preset, include_smooth, scheme, write_vtk, vtk_dir)

Run Riemann cfg6 + reduced DMR (+ optional isentropic vortex order) for one method.
Returns (cases, overall_pass, hard_fails, states_for_vtk).
`states_for_vtk` is Dict name => (state, eq) when write_vtk requested for shock cases.
"""
function run_confirm_suite(
    method_name::AbstractString;
    preset::AbstractString="confirm",
    include_smooth::Bool=true,
    scheme::SchemeConfig=DEFAULT_SCHEME,
    write_vtk::Bool=false,
    vtk_dir::Union{Nothing,AbstractString}=nothing,
)
    spec = get_confirm_preset(preset)
    method = get_capturing_method(method_name)
    cases = Any[]
    hard_fails = String[]
    vtk_states = Dict{String,Any}()

    # --- Isentropic vortex (NullCapturing): multi-D smooth order of the FR core ---
    if include_smooth
        c_v = run_isentropic_vortex_order(;
            p=spec.vortex_p,
            n_list=spec.vortex_n_list,
            t_final=spec.vortex_t,
            cfl=spec.vortex_cfl,
            method=NullCapturing(),
            method_name="null",
        )
        c_v["metrics"]["confirm_role"] = "smooth_order_core"
        c_v["metrics"]["confirm_preset"] = String(preset)
        push!(cases, c_v)
        if !get(c_v, "pass", false)
            push!(hard_fails, "isentropic vortex order failed: $(c_v["observed_orders"])")
        end
        if get(c_v, "diverged", false) || get(c_v, "nan_detected", false)
            push!(hard_fails, "isentropic vortex diverged/nan")
        end
    end

    # --- 2D Riemann cfg6 ---
    c_r, st_r, eq_r = run_euler2d_riemann(;
        p=spec.riemann_p,
        nx=spec.riemann_n,
        ny=spec.riemann_n,
        t_final=spec.riemann_t,
        cfl=spec.riemann_cfl,
        config=:cfg6,
        method=method,
        method_name=method_name,
        require_positivity=true,
    )
    c_r["metrics"]["confirm_role"] = "riemann_cfg6"
    c_r["metrics"]["confirm_preset"] = String(preset)
    c_r["metrics"]["scheme"] = scheme_dict(scheme)
    push!(cases, c_r)
    if !get(c_r, "pass", false)
        push!(hard_fails, "riemann cfg6 failed")
    end
    if get(c_r, "diverged", false) || get(c_r, "nan_detected", false)
        push!(hard_fails, "riemann cfg6 diverged/nan")
    end
    if write_vtk
        vtk_states["riemann_cfg6"] = (st_r, eq_r, method)
    end

    # --- Reduced Double Mach ---
    c_d, st_d, eq_d = run_double_mach_reflection(;
        p=spec.dmr_p,
        nx=spec.dmr_nx,
        ny=spec.dmr_ny,
        t_final=spec.dmr_t,
        cfl=spec.dmr_cfl,
        Lx=spec.dmr_Lx,
        Ly=spec.dmr_Ly,
        strength=:reduced,
        method=method,
        method_name=method_name,
        require_positivity=false,
    )
    c_d["metrics"]["confirm_role"] = "double_mach_reduced"
    c_d["metrics"]["confirm_preset"] = String(preset)
    c_d["metrics"]["scheme"] = scheme_dict(scheme)
    push!(cases, c_d)
    if !get(c_d, "pass", false)
        push!(hard_fails, "double_mach reduced failed")
    end
    if get(c_d, "diverged", false) || get(c_d, "nan_detected", false)
        push!(hard_fails, "double_mach reduced diverged/nan")
    end
    if write_vtk
        vtk_states["double_mach"] = (st_d, eq_d, method)
    end

    overall = isempty(hard_fails) && all(c -> get(c, "pass", false) === true, cases)

    if write_vtk && vtk_dir !== nothing
        mkpath(vtk_dir)
        for (tag, (st, eq, mth)) in pairs(vtk_states)
            path = joinpath(vtk_dir, "$(tag)_$(method_name).vtu")
            write_vtu_high_order_with_capturing(path, st, eq, mth)
            println("  VTU → $path")
        end
    end

    return cases, overall, hard_fails, vtk_states
end

"""
    run_confirm_report(method_name; preset, include_smooth, scheme, write_vtk, vtk_dir)

Build a schema v1 report for one method's confirm suite.
"""
function run_confirm_report(
    method_name::AbstractString;
    preset::AbstractString="confirm",
    include_smooth::Bool=true,
    scheme::SchemeConfig=DEFAULT_SCHEME,
    write_vtk::Bool=false,
    vtk_dir::Union{Nothing,AbstractString}=nothing,
    baseline_name=nothing,
)
    t0 = time()
    cases, overall, hard_fails, _ = run_confirm_suite(
        method_name;
        preset=preset,
        include_smooth=include_smooth,
        scheme=scheme,
        write_vtk=write_vtk,
        vtk_dir=vtk_dir,
    )
    diverged = any(c -> get(c, "diverged", false) === true, cases)
    nan_detected = any(c -> get(c, "nan_detected", false) === true, cases)
    m = get_capturing_method(method_name)
    report = report_skeleton(;
        command="confirm",
        suite="confirm_$(preset)",
        method_name=String(method_name),
        method_params=method_params(m),
        baseline_name=baseline_name,
        overall_pass=overall,
        diverged=diverged,
        nan_detected=nan_detected,
        wall_time_sec=time() - t0,
        hard_gate_failures=hard_fails,
        cases=cases,
        fill_scores=true,
        scheme=scheme,
    )
    report["confirm_preset"] = String(preset)
    report["mesh_summary"] = mesh_summary_dict(
        get_confirm_preset(preset), preset; include_smooth=include_smooth,
    )
    report["scoring_note"] =
        "Confirm scores are diagnostic on fine multi-D meshes; invent composite (coarse 1D quant) remains the historical score of record."
    return report
end

"""True if any case (or the report summary) reports divergence / NaN."""
function _report_diverged_or_nan(report::AbstractDict, cases)
    return get(report, "diverged", false) === true ||
           get(report, "nan_detected", false) === true ||
           any(c -> get(c, "diverged", false) || get(c, "nan_detected", false), cases)
end

"""Pass flag for the first case whose metrics.confirm_role matches `role`."""
function _confirm_case_role_pass(cases, role::AbstractString)
    for c in cases
        met = get(c, "metrics", Dict())
        if string(get(met, "confirm_role", "")) == role
            return get(c, "pass", false) === true
        end
    end
    return false
end

function _confirm_rule_baseline_finished(baseline_report, b_cases)
    bas_div = _report_diverged_or_nan(baseline_report, b_cases)
    bas_overall = get(baseline_report, "overall_pass", false) === true
    bas_ok = !bas_div && bas_overall
    return bas_ok, bas_overall, Dict(
        "id" => "baseline_finished",
        "ok" => bas_ok,
        "notes" => bas_ok ? "baseline completed on fine mesh" :
                   "baseline diverged/nan or overall_pass=false — cannot credit method",
    )
end

function _confirm_rule_method_hard_gates(method_report, m_cases)
    met_div = _report_diverged_or_nan(method_report, m_cases)
    met_overall = get(method_report, "overall_pass", false) === true
    met_hard = !met_div && met_overall
    fails = get(method_report, "hard_gate_failures", String[])
    return met_hard, met_overall, Dict(
        "id" => "method_hard_gates",
        "ok" => met_hard,
        "notes" => met_hard ? "method hard gates passed" :
                   "method failures: $(join(fails, "; "))",
    )
end

function _confirm_rule_multid(m_cases, b_cases)
    r_m = _confirm_case_role_pass(m_cases, "riemann_cfg6")
    r_b = _confirm_case_role_pass(b_cases, "riemann_cfg6")
    d_m = _confirm_case_role_pass(m_cases, "double_mach_reduced")
    d_b = _confirm_case_role_pass(b_cases, "double_mach_reduced")
    multi_ok = r_m && r_b && d_m && d_b
    return Dict(
        "id" => "multid_both_pass",
        "ok" => multi_ok,
        "notes" => "riemann method=$r_m baseline=$r_b; dmr method=$d_m baseline=$d_b",
    )
end

function _confirm_rule_order(m_cases, b_cases, m_scores, b_scores)
    has_smooth =
        any(c -> get(c, "case_type", "") == "smooth_order", m_cases) &&
        any(c -> get(c, "case_type", "") == "smooth_order", b_cases)
    if !has_smooth
        return Dict(
            "id" => "order_no_regression",
            "ok" => true,
            "notes" => "smooth order skipped (--no-smooth)",
        )
    end
    So_m = float(get(m_scores, "order_preservation", 0.0))
    So_b = float(get(b_scores, "order_preservation", 0.0))
    m_order_cases = filter(c -> get(c, "case_type", "") == "smooth_order", m_cases)
    order_pass_m = all(c -> get(c, "order_pass", false) === true, m_order_cases)
    order_ok = (So_m >= So_b - 1e-12) && order_pass_m
    return Dict(
        "id" => "order_no_regression",
        "ok" => order_ok,
        "notes" => "S_order method=$So_m baseline=$So_b order_pass=$order_pass_m",
    )
end

function _confirm_rule_robustness(m_scores, b_scores)
    Sr_m = float(get(m_scores, "robustness", 1.0))
    Sr_b = float(get(b_scores, "robustness", 1.0))
    rob_ok = Sr_m + 1e-12 >= Sr_b
    return Dict(
        "id" => "robustness_no_regression",
        "ok" => rob_ok,
        "notes" => "S_robust method=$Sr_m baseline=$Sr_b",
    )
end

"""
    classify_confirm(method_report, baseline_report; preset) -> Dict

Hard gates + competitiveness on the same fine mesh.
Baseline divergence/NaN ⇒ confirmation_failed (no accidental method "win").
"""
function classify_confirm(
    method_report::AbstractDict,
    baseline_report::AbstractDict;
    preset::AbstractString="confirm",
)
    m_cases = get(method_report, "cases", Any[])
    b_cases = get(baseline_report, "cases", Any[])
    m_scores = get(get(method_report, "summary", Dict()), "scores", Dict())
    b_scores = get(get(baseline_report, "summary", Dict()), "scores", Dict())

    _, bas_overall, r_bas = _confirm_rule_baseline_finished(baseline_report, b_cases)
    _, met_overall, r_met = _confirm_rule_method_hard_gates(method_report, m_cases)
    rules = Dict{String,Any}[
        r_bas,
        r_met,
        _confirm_rule_multid(m_cases, b_cases),
        _confirm_rule_order(m_cases, b_cases, m_scores, b_scores),
        _confirm_rule_robustness(m_scores, b_scores),
    ]

    all_ok = all(r -> r["ok"] === true, rules)
    status = all_ok ? "confirmed" : "confirmation_failed"

    return Dict{String,Any}(
        "confirmation_status" => status,
        "preset" => String(preset),
        "baseline_name" => get(baseline_report, "method_name", nothing),
        "method_name" => get(method_report, "method_name", nothing),
        "method_overall_pass" => met_overall,
        "baseline_overall_pass" => bas_overall,
        "rules" => rules,
        "absolute_scores" => m_scores,
        "baseline_absolute_scores" => b_scores,
        "mesh_summary" => get(method_report, "mesh_summary", nothing),
        "scoring_note" => get(method_report, "scoring_note", nothing),
    )
end

function print_confirm_summary(method_name::AbstractString, cmp::AbstractDict; io::IO=stdout)
    status = cmp["confirmation_status"]
    println(io, "Confirm: $method_name    status: $status    preset: $(get(cmp, "preset", "?"))")
    println(
        io,
        "  method_overall_pass=$(cmp["method_overall_pass"])  baseline_overall_pass=$(cmp["baseline_overall_pass"])",
    )
    for r in get(cmp, "rules", [])
        mark = get(r, "ok", false) ? "OK" : "FAIL"
        println(io, "  [$mark] $(r["id"]): $(r["notes"])")
    end
    return nothing
end

"""
    entry_from_confirm(method_name, method_report, baseline_report, cmp; ...) -> Dict

Build experiment-log entry for a confirm run.
"""
function entry_from_confirm(
    method_name::AbstractString,
    method_report::AbstractDict,
    baseline_report::AbstractDict,
    cmp::AbstractDict;
    hypothesis::AbstractString="",
    lessons::AbstractString="",
    strengths::AbstractString="",
    weaknesses::AbstractString="",
    scheme=FROZEN_INVENT_SCHEME,
    git_ref::AbstractString="",
    artifacts::Union{Nothing,AbstractDict}=nothing,
    entry_id::Union{Nothing,AbstractString}=nothing,
    date::Date=Dates.today(),
    preset::AbstractString="confirm",
)
    status = String(get(cmp, "confirmation_status", "confirmation_failed"))
    scores = get(cmp, "absolute_scores", Dict())
    mesh = get(cmp, "mesh_summary", get(method_report, "mesh_summary", Dict()))
    mesh_note = if mesh isa AbstractDict && haskey(mesh, "preset")
        mesh_note_string(get_confirm_preset(string(mesh["preset"])), string(mesh["preset"]);
            include_smooth=haskey(mesh, "vortex"))
    else
        "preset=$preset"
    end

    hyp = String(hypothesis)
    les = String(lessons)
    isempty(strip(hyp)) &&
        (hyp = "Fine-mesh confirmation of short-listed method (coarse invent remains score of record).")
    isempty(strip(les)) &&
        (les = status == "confirmed" ?
         "Confirmed on fine multi-D meshes; safe to consider robustness + snapshot freeze with --require-confirm." :
         "Confirmation failed — do not promote; inspect rules/artifacts before retrying.")

    arts = if artifacts !== nothing
        Dict{String,Any}(String(k) => v for (k, v) in pairs(artifacts))
    else
        Dict{String,Any}()
    end

    sch = if scheme isa NamedTuple
        Dict{String,Any}(
            "points" => string(scheme.points),
            "flux" => string(scheme.flux),
            "time" => string(scheme.time),
        )
    elseif scheme isa SchemeConfig
        scheme_dict(scheme)
    else
        Dict{String,Any}(
            "points" => string(get(scheme, "points", FROZEN_INVENT_SCHEME.points)),
            "flux" => string(get(scheme, "flux", FROZEN_INVENT_SCHEME.flux)),
            "time" => string(get(scheme, "time", FROZEN_INVENT_SCHEME.time)),
        )
    end

    suffix = preset == "confirm" ? "confirm" : "confirm_$(preset)"
    return Dict{String,Any}(
        "id" => something(entry_id, make_entry_id(method_name; date=date, suffix=suffix)),
        "date" => string(date),
        "method" => String(method_name),
        "baseline" => string(something(get(cmp, "baseline_name", nothing), "persson_av")),
        "hypothesis" => hyp,
        "scheme" => sch,
        "metrics" => Dict{String,Any}(
            "confirmation_status" => status,
            "confirm_preset" => String(preset),
            "mesh" => mesh_note,
            "order_preservation" => get(scores, "order_preservation", nothing),
            "dissipation" => get(scores, "dissipation", nothing),
            "shock_quality" => get(scores, "shock_quality", nothing),
            "robustness" => get(scores, "robustness", nothing),
            "composite" => get(scores, "composite", nothing),
            "method_overall_pass" => get(cmp, "method_overall_pass", nothing),
            "baseline_overall_pass" => get(cmp, "baseline_overall_pass", nothing),
        ),
        "strengths" => String(strengths),
        "weaknesses" => String(weaknesses),
        "lessons" => les,
        "status" => status,
        "artifacts" => arts,
        "git_ref" => String(git_ref),
        "narrative_complete" => true,
        "kind" => "confirm",
    )
end

"""
    confirm_method(method_name; baseline, preset, ...) -> (met, bas, cmp)

Orchestrate fine-mesh confirmation for method vs baseline.
"""
function confirm_method(
    method_name::AbstractString;
    baseline::AbstractString="persson_av",
    preset::AbstractString="confirm",
    report_dir::AbstractString="results/confirm",
    include_smooth::Bool=true,
    write_vtk::Bool=false,
    scheme::SchemeConfig=DEFAULT_SCHEME,
    append_log::Bool=true,
    log_path::Union{Nothing,AbstractString}=nothing,
    yaml_path::Union{Nothing,AbstractString}=nothing,
    hypothesis::AbstractString="",
    lessons::AbstractString="",
    strengths::AbstractString="",
    weaknesses::AbstractString="",
    git_ref::AbstractString="",
    threads::Int=1,
)
    method_name = require_registered_method(method_name)
    baseline = require_registered_method(baseline; role="baseline")
    spec = get_confirm_preset(preset)
    mkpath(report_dir)
    vtk_dir = write_vtk ? joinpath(report_dir, "vtu") : nothing

    # Official confirm-for-promotion is always serial residual.
    # threads>1 is exploratory only and does NOT satisfy publication_grade confirm.
    use_threads = max(1, Int(threads))
    official = use_threads <= 1

    println("=== FRForge confirm (fine-mesh short-list gate) ===")
    println("method=$method_name  baseline=$baseline  preset=$preset")
    println("scheme=", scheme_dict(scheme))
    println("mesh: ", mesh_note_string(spec, preset; include_smooth=include_smooth))
    println(spec.wall_time_note)
    println(
        "Note: invent composite history is unchanged; this is a separate confirmation stream.",
    )
    if !official
        println()
        println(
            "WARNING: --threads=$use_threads requested. Threaded confirm is INFORMATIONAL ONLY.",
        )
        println(
            "  It does NOT satisfy fine-mesh confirmation for publication_grade or",
        )
        println(
            "  `snapshot freeze --require-confirm`. Re-run with serial residual (default) for promotion.",
        )
        println()
    end
    if scheme != DEFAULT_SCHEME
        println(
            "WARNING: non-default scheme — confirm result does not count for publication_grade.",
        )
    end

    return with_frforge_threads(use_threads) do
        paths = report_trio_paths(report_dir, method_name, baseline; tag=String(preset))
        bas_path, met_path, cmp_path = paths.baseline, paths.method, paths.compare

        println("residual threads: $(frforge_thread_count())  (julia nthreads=$(Threads.nthreads()))")
        println("\n--- Baseline confirm: $baseline ---")
        bas = run_confirm_report(
            baseline;
            preset=preset,
            include_smooth=include_smooth,
            scheme=scheme,
            write_vtk=false,
            baseline_name=nothing,
        )
        stamp_workflow_report!(
            bas;
            command="confirm",
            baseline_name=nothing,
            scheme=scheme,
            extra=Dict{String,Any}("confirm_official" => official),
        )
        write_report(bas_path, bas)

        println("\n--- Method confirm: $method_name ---")
        met = run_confirm_report(
            method_name;
            preset=preset,
            include_smooth=include_smooth,
            scheme=scheme,
            write_vtk=write_vtk,
            vtk_dir=vtk_dir,
            baseline_name=baseline,
        )
        stamp_workflow_report!(
            met;
            command="confirm",
            baseline_name=baseline,
            scheme=scheme,
            extra=Dict{String,Any}("confirm_official" => official),
        )
        write_report(met_path, met)

        cmp = classify_confirm(met, bas; preset=preset)
        if !official
            # Do not claim confirmed for promotion when threaded
            cmp["confirmation_status"] = "confirmation_failed"
            cmp["informational_threaded"] = true
            push!(
                cmp["rules"],
                Dict(
                    "id" => "serial_residual_required",
                    "ok" => false,
                    "notes" => "threaded confirm is informational only; not valid for publication_grade",
                ),
            )
        end
        met["confirmation_status"] = cmp["confirmation_status"]
        write_report(met_path, met)
        write_json_pretty(cmp_path, cmp)

        print_confirm_summary(method_name, cmp)
        print_report_trio(met_path, bas_path, cmp_path)

        if append_log
            lp = something(log_path, default_experiment_log_path())
            yp = yaml_path === nothing ? default_experiment_log_yaml_path() : yaml_path
            yp_use = (yp isa AbstractString && isempty(yp)) ? nothing : yp
            arts = Dict{String,Any}(
                "method_report" => met_path,
                "baseline_report" => bas_path,
                "compare" => cmp_path,
            )
            if write_vtk && vtk_dir !== nothing
                arts["vtu_dir"] = vtk_dir
            end
            entry = entry_from_confirm(
                method_name,
                met,
                bas,
                cmp;
                hypothesis=hypothesis,
                lessons=lessons,
                strengths=strengths,
                weaknesses=weaknesses,
                scheme=scheme,
                git_ref=git_ref,
                artifacts=arts,
                preset=preset,
            )
            if !official
                entry["status"] = "confirmation_failed"
                entry["metrics"]["confirmation_status"] = "confirmation_failed"
                entry["metrics"]["informational_threaded"] = true
                entry["lessons"] =
                    "Threaded confirm only — re-run serial confirm before promotion. " *
                    string(get(entry, "lessons", ""))
            end
            append_experiment_entry!(entry; path=lp, yaml_path=yp_use)
            println("Experiment log appended: $(entry["id"])  →  $lp")
        end

        return met, bas, cmp
    end
end

"""
    method_has_confirm_pass(method; log_path, confirm_compare) -> (Bool, String)

True if a confirmed log entry or confirm compare JSON exists for method.
"""
function method_has_confirm_pass(
    method::AbstractString;
    log_path::AbstractString=default_experiment_log_path(),
    confirm_compare::Union{Nothing,AbstractString}=nothing,
)
    if confirm_compare !== nothing && isfile(confirm_compare)
        cmp = JSON.parsefile(confirm_compare)
        if string(get(cmp, "confirmation_status", "")) == "confirmed"
            return true, "confirm compare JSON: $confirm_compare"
        end
        return false, "confirm compare present but not confirmed: $confirm_compare"
    end
    if isfile(log_path)
        entries = parse_experiment_log(log_path)
        for e in reverse(entries)
            String(get(e, "method", "")) != String(method) && continue
            st = String(get(e, "status", ""))
            met = get(e, "metrics", Dict())
            cs = string(get(met, "confirmation_status", ""))
            if st == "confirmed" || cs == "confirmed"
                return true, "log entry $(e["id"]) status=confirmed"
            end
        end
    end
    return false, "no confirmed log entry or compare JSON for method=$method"
end
