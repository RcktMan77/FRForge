# `frforge snapshot` freeze/verify/tables.

function _parse_snapshot_args(args)
    if isempty(args) || args[1] in ("-h", "--help", "help")
        return Dict("sub" => "help")
    end
    sub = args[1]
    rest = args[2:end]
    if sub == "freeze"
        s = ArgParseSettings(description="Freeze a reproducibility snapshot.", prog="frforge snapshot freeze")
        @add_arg_table! s begin
            "--method", "-m"
            required = true
            "--baseline", "-b"
            default = "persson_av"
            "--method-report"
            required = true
            dest_name = "method_report"
            "--baseline-report"
            required = true
            dest_name = "baseline_report"
            "--compare"
            default = ""
            "--out"
            help = "Snapshot root directory"
            default = "results/snapshots"
            "--git-ref"
            default = ""
            dest_name = "git_ref"
            "--source"
            help = "Extra source file (repeatable escape hatch)"
            action = :append_arg
            default = String[]
            "--no-append-log"
            action = :store_true
            dest_name = "no_append_log"
            "--log-path"
            default = ""
            dest_name = "log_path"
            "--require-confirm"
            help = "Hard-fail freeze unless fine-mesh confirm passed (recommended for papers)"
            action = :store_true
            dest_name = "require_confirm"
            "--confirm-compare"
            help = "Optional path to confirm compare JSON"
            default = ""
            dest_name = "confirm_compare"
        end
        opts = parse_args(rest, s)
        opts["sub"] = "freeze"
        return opts
    elseif sub == "verify"
        s = ArgParseSettings(
            description="Verify snapshot (default: cheap hash check; --rerun for invent).",
            prog="frforge snapshot verify",
        )
        @add_arg_table! s begin
            "path"
            help = "Snapshot directory"
            required = true
            "--rerun"
            help = "Explicit full invent re-run (slow; not for required CI)"
            action = :store_true
            "--require-git-ref"
            action = :store_true
            dest_name = "require_git_ref"
            "--tol-rel"
            arg_type = Float64
            default = 1e-10
            dest_name = "tol_rel"
            "--tol-abs"
            arg_type = Float64
            default = 1e-12
            dest_name = "tol_abs"
        end
        opts = parse_args(rest, s)
        opts["sub"] = "verify"
        return opts
    elseif sub == "tables"
        s = ArgParseSettings(description="Tables from frozen snapshot JSON.", prog="frforge snapshot tables")
        @add_arg_table! s begin
            "path"
            required = true
            "--out"
            help = "Markdown output path"
            default = ""
            "--csv"
            help = "CSV output path"
            default = ""
        end
        opts = parse_args(rest, s)
        opts["sub"] = "tables"
        return opts
    else
        return Dict("sub" => "unknown", "cmd" => sub)
    end
end

function cli_snapshot(args)
    opts = _parse_snapshot_args(args)
    sub = opts["sub"]
    if sub == "help"
        println("Usage: frforge snapshot {freeze|verify|tables}")
        println("  freeze   Package method sources + reports + hashes (immutable)")
        println("  verify   Cheap hash/manifest check (default); add --rerun for invent")
        println("  tables   Markdown/CSV comparison tables from frozen JSON")
        println("Freeze only after invent short-list + fine-mesh confirm (+ robustness).")
        println("Paper-facing: add --require-confirm after `frforge confirm` succeeds.")
        return 0
    elseif sub == "freeze"
        compare = isempty(opts["compare"]) ? nothing : opts["compare"]
        log_path = isempty(opts["log_path"]) ? default_experiment_log_path() : opts["log_path"]
        conf_cmp = isempty(opts["confirm_compare"]) ? nothing : opts["confirm_compare"]
        extra = String.(opts["source"])
        dir = freeze_snapshot(;
            method=opts["method"],
            baseline=opts["baseline"],
            method_report=opts["method_report"],
            baseline_report=opts["baseline_report"],
            compare=compare,
            out_root=opts["out"],
            git_ref=opts["git_ref"],
            source_extra=extra,
            log_path=log_path,
            append_log=!opts["no_append_log"],
            require_confirm=opts["require_confirm"],
            confirm_compare=conf_cmp,
        )
        println("Snapshot frozen: $dir")
        return 0
    elseif sub == "verify"
        res = verify_snapshot(
            opts["path"];
            rerun=opts["rerun"],
            tol_rel=opts["tol_rel"],
            tol_abs=opts["tol_abs"],
            require_git_ref=opts["require_git_ref"],
        )
        println("mode=$(res["mode"]) ok=$(res["ok"]) method=$(res["method"])")
        if !isempty(res["errors"])
            for e in res["errors"]
                println(stderr, "  ERROR: ", e)
            end
        end
        return res["ok"] ? 0 : 1
    elseif sub == "tables"
        out_md = isempty(opts["out"]) ? nothing : opts["out"]
        out_csv = isempty(opts["csv"]) ? nothing : opts["csv"]
        tab = snapshot_tables(opts["path"]; out_md=out_md, out_csv=out_csv)
        print(tab["markdown"])
        out_md !== nothing && println("Wrote $out_md")
        out_csv !== nothing && println("Wrote $out_csv")
        return 0
    else
        println(stderr, "Unknown snapshot subcommand: ", get(opts, "cmd", sub))
        println(stderr, "  Use: frforge snapshot {$(join(CLI_SNAPSHOT_SUBCOMMANDS, "|"))}")
        return 2
    end
end

