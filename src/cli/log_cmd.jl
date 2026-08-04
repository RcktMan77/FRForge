# `frforge log` subcommands.

function _parse_log_args(args)
    if isempty(args) || args[1] in ("-h", "--help", "help")
        return Dict("sub" => "help")
    end
    sub = args[1]
    rest = args[2:end]
    if sub == "list"
        s = ArgParseSettings(description="List experiment log entry ids.", prog="frforge log list")
        @add_arg_table! s begin
            "--path"
            help = "Path to experiment_log.md"
            default = ""
        end
        opts = parse_args(rest, s)
        opts["sub"] = "list"
        return opts
    elseif sub == "summary"
        s = ArgParseSettings(description="Summarize experiment log.", prog="frforge log summary")
        @add_arg_table! s begin
            "--path"
            help = "Path to experiment_log.md"
            default = ""
            "--json"
            help = "Optional path to write machine-readable JSON"
            default = ""
        end
        opts = parse_args(rest, s)
        opts["sub"] = "summary"
        return opts
    elseif sub == "frontier"
        s = ArgParseSettings(description="Frontier / near-miss methods.", prog="frforge log frontier")
        @add_arg_table! s begin
            "--path"
            default = ""
            "--json"
            default = ""
        end
        opts = parse_args(rest, s)
        opts["sub"] = "frontier"
        return opts
    elseif sub == "lessons"
        s = ArgParseSettings(description="Index of lessons and weaknesses.", prog="frforge log lessons")
        @add_arg_table! s begin
            "--path"
            default = ""
            "--query", "-q"
            help = "Case-insensitive substring filter"
            default = ""
            "--json"
            default = ""
        end
        opts = parse_args(rest, s)
        opts["sub"] = "lessons"
        return opts
    elseif sub == "show"
        s = ArgParseSettings(description="Show one log entry by id.", prog="frforge log show")
        @add_arg_table! s begin
            "id"
            help = "Entry id (### heading)"
            required = true
            "--path"
            default = ""
            "--json"
            default = ""
        end
        opts = parse_args(rest, s)
        opts["sub"] = "show"
        return opts
    elseif sub == "pareto"
        s = ArgParseSettings(description="Pareto-style order/dissip/shock table.", prog="frforge log pareto")
        @add_arg_table! s begin
            "--path"
            default = ""
            "--json"
            default = ""
        end
        opts = parse_args(rest, s)
        opts["sub"] = "pareto"
        return opts
    elseif sub == "append"
        s = ArgParseSettings(
            description="Append experiment log entry from two invent JSON reports.",
            prog="frforge log append",
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
            "--path"
            help = "Path to experiment_log.md"
            default = ""
            "--hypothesis"
            default = ""
            "--lessons"
            default = ""
            "--delta"
            arg_type = Float64
            default = DEFAULT_SCORE_MARGIN
            "--vtk-produced"
            action = :store_true
            dest_name = "vtk_produced"
        end
        opts = parse_args(rest, s)
        opts["sub"] = "append"
        return opts
    else
        return Dict("sub" => "unknown", "cmd" => sub)
    end
end

function _log_path(opts)
    p = get(opts, "path", "")
    return isempty(p) ? default_experiment_log_path() : p
end

function _maybe_write_json(path::AbstractString, obj)
    isempty(path) && return nothing
    write_json_pretty(path, obj)
    println("Wrote JSON → $path")
    return nothing
end

function cli_log(args)
    opts = _parse_log_args(args)
    sub = opts["sub"]
    if sub == "help"
        println("Usage: frforge log {$(join(CLI_LOG_SUBCOMMANDS, "|"))} [options]")
        println("  list      List entry ids")
        println("  summary   Counts, latest per method, narrative TODOs")
        println("  frontier  Baseline + promising-class + near-miss pass_gates")
        println("  pareto    Order/dissipation/shock table with non-dominated flags")
        println("  lessons   Flatten lessons/weaknesses (--query filter)")
        println("  show <id> Print one full entry")
        println("  append    Append from invent JSON reports")
        println("Common: --path log.md  --json out.json")
        return 0
    elseif sub == "list"
        path = _log_path(opts)
        ids = list_experiment_entry_ids(path)
        println("Experiment log: $path")
        println("n_entries=$(length(ids))")
        for id in ids
            println("  - ", id)
        end
        return 0
    elseif sub == "summary"
        entries = parse_experiment_log(_log_path(opts))
        summary = log_summary(entries)
        print(format_log_summary_text(summary))
        _maybe_write_json(get(opts, "json", ""), summary)
        return 0
    elseif sub == "frontier"
        entries = parse_experiment_log(_log_path(opts))
        rows = log_frontier(entries)
        print(format_log_frontier_text(rows))
        _maybe_write_json(get(opts, "json", ""), rows)
        return 0
    elseif sub == "pareto"
        entries = parse_experiment_log(_log_path(opts))
        rows = log_pareto(entries)
        print(format_log_pareto_text(rows))
        _maybe_write_json(get(opts, "json", ""), rows)
        return 0
    elseif sub == "lessons"
        entries = parse_experiment_log(_log_path(opts))
        q = get(opts, "query", "")
        rows = log_lessons(entries; query=isempty(q) ? nothing : q)
        print(format_log_lessons_text(rows))
        _maybe_write_json(get(opts, "json", ""), rows)
        return 0
    elseif sub == "show"
        entries = parse_experiment_log(_log_path(opts))
        e = get_experiment_entry(entries, opts["id"])
        if e === nothing
            println(stderr, "Entry not found: ", opts["id"])
            return 1
        end
        print(format_log_entry_text(e))
        _maybe_write_json(get(opts, "json", ""), e)
        return 0
    elseif sub == "append"
        path = _log_path(opts)
        met = load_report(opts["method_report"])
        bas = load_report(opts["baseline_report"])
        cmp = classify_candidate(met, bas; δ=opts["delta"], vtk_produced=opts["vtk_produced"])
        method_name = string(get(met, "method_name", "method"))
        arts = Dict{String,Any}(
            "method_report" => opts["method_report"],
            "baseline_report" => opts["baseline_report"],
        )
        entry = invent_append_log!(
            method_name,
            met,
            bas,
            cmp;
            log_path=path,
            yaml_path=nothing,
            artifacts=arts,
            hypothesis=opts["hypothesis"],
            lessons=opts["lessons"],
        )
        println("Appended $(entry["id"]) → $path")
        return 0
    else
        println(stderr, "Unknown log subcommand: ", get(opts, "cmd", sub))
        println(stderr, "  Use: frforge log {$(join(CLI_LOG_SUBCOMMANDS, "|"))}")
        return 2
    end
end

