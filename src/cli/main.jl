# CLI process entry: `main_cli` dispatcher.
#
# Command implementations live in sibling files (see constants.jl layout note).

"""
    main_cli(args=ARGS) -> Int

Top-level CLI dispatcher. Returns a process exit code.
"""
function main_cli(args = ARGS)
    if isempty(args) || args[1] in ("-h", "--help", "help")
        _print_usage()
        return isempty(args) ? 2 : 0
    end

    cmd = args[1]
    rest = args[2:end]

    try
        if cmd == "test"
            return cli_test(_parse_test_args(rest))
        elseif cmd == "run"
            return cli_run(rest)
        elseif cmd == "invent"
            return cli_invent(rest)
        elseif cmd == "confirm"
            return cli_confirm(rest)
        elseif cmd == "score"
            return cli_score(rest)
        elseif cmd == "log"
            return cli_log(rest)
        elseif cmd == "snapshot"
            return cli_snapshot(rest)
        elseif cmd == "robustness"
            return cli_robustness(rest)
        else
            println(stderr, "Unknown command: $cmd")
            println(stderr, "  Valid commands: $(join(CLI_COMMANDS, ", "))")
            _print_usage()
            return 2
        end
    catch e
        if e isa ArgParseError
            println(stderr, "Argument error: ", e.text)
            return 2
        end
        rethrow()
    end
end
