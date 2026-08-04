# CLI shared constants and top-level usage text.
#
# Layout (cli/):
#   scheme_args.jl  — --points/--flux/--time helpers
#   constants.jl    — command/suite/case lists + usage
#   run.jl          — frforge run (per-case helpers)
#   test_cmd.jl     — frforge test
#   invent_cmd.jl   — invent / score / confirm / robustness
#   log_cmd.jl      — experiment log subcommands
#   snapshot_cmd.jl — snapshot freeze/verify/tables
#   main.jl         — main_cli dispatcher

const CLI_COMMANDS = ("test", "run", "invent", "score", "confirm", "log", "snapshot", "robustness")

# Keep help strings and error messages in lockstep with `_dispatch_test_suite`.
const CLI_TEST_SUITE_HELP =
    "Suite: smoke|advection|burgers|euler|capturing|quant|2d|2d_capturing|curved|benchmarks|optional2d|full " *
    "(aliases: m1–m5, m8, p31–p33b, riemann, vortex, dmr, ffs)"

const CLI_RUN_CASES = (
    "advection_sine",
    "burgers_square",
    "euler_density_wave",
    "sod",
    "shu_osher",
    "advection2d",
    "euler2d_wave",
    "euler2d_jump",
    "riemann2d",
    "double_mach",
)

const CLI_RUN_CASE_HELP = "Case: " * join(CLI_RUN_CASES, "|")

const CLI_LOG_SUBCOMMANDS = ("list", "summary", "frontier", "lessons", "show", "pareto", "append")
const CLI_SNAPSHOT_SUBCOMMANDS = ("freeze", "verify", "tables")

function _print_usage(io=stderr)
    println(io, "FRForge — high-order Flux Reconstruction laboratory")
    println(io, "Usage: frforge {$(join(CLI_COMMANDS, "|"))} [options]")
    println(io)
    println(io, "Commands:")
    println(io, "  test         Run verification suite and emit JSON report")
    println(io, "  run          Run a single case")
    println(io, "  invent       Coarse quant suite + classify (fast; score history)")
    println(io, "  confirm      Fine-mesh multi-D confirmation after short-list")
    println(io, "  score        Classify two existing JSON reports (method vs baseline)")
    println(io, "  log          Experiment log: $(join(CLI_LOG_SUBCOMMANDS, "|"))")
    println(io, "  snapshot     Freeze/verify/tables for reproducibility packages")
    println(io, "  robustness   Scheme robustness matrix (CI-light or full/nightly)")
    println(io)
    println(io, "Agent workflow: invent (coarse) → log frontier → confirm (fine) → snapshot freeze.")
    println(
        io,
        "Frozen invent scheme: $(FROZEN_INVENT_SCHEME.points) + $(FROZEN_INVENT_SCHEME.flux) + $(FROZEN_INVENT_SCHEME.time)",
    )
    println(io, "Official invent/confirm/scoring use serial residual (bit-deterministic).")
    println(io)
    println(io, "Examples:")
    println(io, "  frforge test --suite quant --method persson_av --report results/quant/report.json")
    println(io, "  frforge invent --method scaled_persson --baseline persson_av")
    println(io, "  frforge confirm --method scaled_persson --baseline persson_av")
    println(io, "  frforge confirm --method scaled_persson --preset quick   # smoke")
    println(io, "  frforge score --method-report a.json --baseline-report b.json")
    println(io, "  frforge log list")
    println(io, "  frforge log summary")
    println(io, "  frforge log frontier")
    println(io, "  frforge log lessons --query c_av")
    println(io, "  frforge robustness --method scaled_persson --matrix ci")
    println(io, "  frforge snapshot freeze ... --require-confirm   # paper-facing")
    println(io, "  frforge run --case sod --p 2 --ne 64 --method persson_av")
    println(io, describe_methods())
end

