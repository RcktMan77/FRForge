# `frforge invent|score|confirm|robustness` CLI entry points.

function _parse_invent_args(args)
    s = ArgParseSettings(
        description = "Invent: run quant suite for method vs baseline and classify.",
        prog = "frforge invent",
    )
    @add_arg_table! s begin
        "--method", "-m"
        help = "Candidate method name (registered)"
        required = true
        "--baseline", "-b"
        help = "Baseline method name"
        default = "persson_av"
        "--report-dir"
        help = "Directory for invent JSON artifacts"
        default = "results/invent"
        dest_name = "report_dir"
        "--delta"
        help = "Composite margin δ_score for promising status"
        arg_type = Float64
        default = DEFAULT_SCORE_MARGIN
        "--vtk-produced"
        help = "Set true if HO VTK was produced (accepted_candidate)"
        action = :store_true
        dest_name = "vtk_produced"
        "--no-append-log"
        help = "Do not append an entry to research/experiment_log.md"
        action = :store_true
        dest_name = "no_append_log"
        "--hypothesis"
        help = "Hypothesis text for experiment log entry"
        default = ""
        "--lessons"
        help = "Lessons text for experiment log entry (required for promising+)"
        default = ""
        "--log-path"
        help = "Override experiment log path (default: research/experiment_log.md)"
        default = ""
        dest_name = "log_path"
    end
    return parse_args(args, s)
end

function _parse_score_args(args)
    s = ArgParseSettings(
        description = "Score two existing reports (method vs baseline).",
        prog = "frforge score",
    )
    @add_arg_table! s begin
        "--method-report"
        help = "Path to method report JSON"
        required = true
        dest_name = "method_report"
        "--baseline-report"
        help = "Path to baseline report JSON"
        required = true
        dest_name = "baseline_report"
        "--delta"
        help = "Composite margin δ_score"
        arg_type = Float64
        default = DEFAULT_SCORE_MARGIN
        "--vtk-produced"
        action = :store_true
        dest_name = "vtk_produced"
        "--output", "-o"
        help = "Optional path for comparison JSON"
        default = ""
    end
    return parse_args(args, s)
end

function cli_invent(args)
    opts = _parse_invent_args(args)
    method = require_registered_method(opts["method"])
    baseline = require_registered_method(opts["baseline"]; role = "baseline")
    log_path = isempty(opts["log_path"]) ? nothing : opts["log_path"]
    met, bas, cmp = invent_method(
        method;
        baseline = baseline,
        report_dir = opts["report_dir"],
        δ = opts["delta"],
        vtk_produced = opts["vtk_produced"],
        append_log = !opts["no_append_log"],
        log_path = log_path,
        hypothesis = opts["hypothesis"],
        lessons = opts["lessons"],
    )
    status = cmp["candidate_status"]
    # Exit 0 for promising/accepted/pass_gates; 1 for rejected
    return status == "rejected" ? 1 : 0
end

function cli_score(args)
    opts = _parse_score_args(args)
    out = isempty(opts["output"]) ? nothing : opts["output"]
    cmp = score_reports(
        opts["method_report"],
        opts["baseline_report"];
        δ = opts["delta"],
        vtk_produced = opts["vtk_produced"],
        out_path = out,
    )
    return cmp["candidate_status"] == "rejected" ? 1 : 0
end

function _parse_confirm_args(args)
    s = ArgParseSettings(
        description =
        "Fine-mesh confirmation after invent short-list. " *
        "Default preset=confirm (~10–30 min class: Riemann 64² p=2 + DMR 120×40 + vortex). " *
        "Use --preset presentation for paper meshes (hours). --preset quick for smoke. " *
        "Does not change invent composite-score history.",
        prog = "frforge confirm",
    )
    @add_arg_table! s begin
        "--method", "-m"
        help = "Candidate method name (registered)"
        required = true
        "--baseline", "-b"
        help = "Baseline method name"
        default = "persson_av"
        "--preset"
        help = "Mesh preset: confirm (default) | presentation | quick"
        default = "confirm"
        "--report-dir"
        help = "Directory for confirm JSON / VTU"
        default = "results/confirm"
        dest_name = "report_dir"
        "--no-smooth"
        help = "Skip isentropic vortex order study"
        action = :store_true
        dest_name = "no_smooth"
        "--vtk"
        help = "Write HO VTU for Riemann/DMR under report-dir/vtu/"
        action = :store_true
        "--no-append-log"
        help = "Do not append experiment log entry"
        action = :store_true
        dest_name = "no_append_log"
        "--hypothesis"
        help = "Hypothesis text for experiment log entry"
        default = ""
        "--lessons"
        help = "Lessons text for experiment log entry"
        default = ""
        "--log-path"
        help = "Override experiment log path (default: research/experiment_log.md)"
        default = ""
        dest_name = "log_path"
        "--points"
        help = "Solution points (default gl; non-default not for publication_grade)"
        default = "gl"
        "--flux"
        help = "Numerical flux: rusanov (default) | hllc (non-default not for publication_grade)"
        default = "rusanov"
        "--time"
        help = "Time integrator: ssp_rk3 (default) | ssp_rk2 (non-default not for publication_grade)"
        default = "ssp_rk3"
        "--threads"
        help = "Residual threads (default 1=serial/official). N>1 is informational only — not valid for publication_grade."
        arg_type = Int
        default = 1
    end
    return parse_args(args, s)
end

function cli_confirm(args)
    opts = _parse_confirm_args(args)
    scheme = parse_scheme(; points = opts["points"], flux = opts["flux"], time = opts["time"])
    log_path = isempty(opts["log_path"]) ? nothing : opts["log_path"]
    _, _, cmp = confirm_method(
        opts["method"];
        baseline = opts["baseline"],
        preset = opts["preset"],
        report_dir = opts["report_dir"],
        include_smooth = !opts["no_smooth"],
        write_vtk = opts["vtk"],
        scheme = scheme,
        append_log = !opts["no_append_log"],
        log_path = log_path,
        hypothesis = opts["hypothesis"],
        lessons = opts["lessons"],
        threads = opts["threads"],
    )
    return cmp["confirmation_status"] == "confirmed" ? 0 : 1
end

function _parse_robustness_args(args)
    s = ArgParseSettings(
        description = "Robustness matrix across scheme axes (points × flux × time).",
        prog = "frforge robustness",
    )
    @add_arg_table! s begin
        "--method", "-m"
        help = "Capturing method name"
        required = true
        "--baseline", "-b"
        help = "Baseline method"
        default = "persson_av"
        "--matrix"
        help = "Matrix size: ci (default, required-CI safe) | full (local/nightly)"
        default = "ci"
        "--report-dir"
        help = "Root directory for robustness JSON"
        default = "results/robustness"
        dest_name = "report_dir"
        "--delta"
        help = "Composite margin δ for invent classification per cell"
        arg_type = Float64
        default = DEFAULT_SCORE_MARGIN
        "--full-suite"
        help = "Force full quant suite even for matrix=ci (slow; not for required CI)"
        action = :store_true
        dest_name = "full_suite"
        "--light"
        help = "Force light suite even for matrix=full"
        action = :store_true
        "--no-append-log"
        help = "Do not append experiment log entry"
        action = :store_true
        dest_name = "no_append_log"
        "--narrative-complete"
        help = "Mark hypothesis/lessons as complete for promotion assessment"
        action = :store_true
        dest_name = "narrative_complete"
        "--invent-status"
        help = "Override invent status for default-scheme promotion rule (e.g. promising)"
        default = ""
        dest_name = "invent_status"
        "--hypothesis"
        default = ""
        "--lessons"
        default = ""
        "--log-path"
        default = ""
        dest_name = "log_path"
    end
    return parse_args(args, s)
end

function cli_robustness(args)
    opts = _parse_robustness_args(args)
    method = opts["method"]
    baseline = opts["baseline"]
    method = require_registered_method(method)
    baseline = require_registered_method(baseline; role = "baseline")
    matrix = Symbol(lowercase(opts["matrix"]))
    matrix in (:ci, :full) || error("--matrix must be ci or full (got $(opts["matrix"]))")
    light = if opts["full_suite"]
        false
    elseif opts["light"]
        true
    else
        nothing  # auto: ci→light, full→full
    end
    invent_status = isempty(opts["invent_status"]) ? nothing : opts["invent_status"]
    log_path = isempty(opts["log_path"]) ? nothing : opts["log_path"]
    summary = run_robustness_matrix(
        method;
        baseline = baseline,
        matrix = matrix,
        report_dir = opts["report_dir"],
        light = light,
        δ = opts["delta"],
        append_log = !opts["no_append_log"],
        log_path = log_path,
        narrative_complete = opts["narrative_complete"],
        invent_status = invent_status,
        hypothesis = opts["hypothesis"],
        lessons = opts["lessons"],
    )
    # Exit 0 if all cells ok (not necessarily publication-grade)
    all_ok = all(c -> get(c, "ok", false) === true, summary["cells"])
    return all_ok ? 0 : 1
end

"""
    main_cli(args=ARGS) -> Int

Top-level CLI dispatcher. Returns a process exit code.
"""
