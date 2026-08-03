"""
    FRForge

High-order Flux Reconstruction laboratory for inventing and evaluating
novel shock-capturing methods on the compressible Euler equations.

This package is intentionally green-field: FR operators, residual evaluation,
verification, and I/O are implemented from first principles with clarity and
verifiability prioritized over performance in early milestones.
"""
module FRForge

using ArgParse
using Dates
using JSON

export main_cli
export write_report_skeleton, load_report, validate_report_keys
export DEFAULT_SCORING_WEIGHTS, SCORING_FORMULA_VERSION, SCHEMA_VERSION
export REQUIRED_TOP_LEVEL_KEYS, REQUIRED_SUMMARY_KEYS

include("verification/schema_keys.jl")
include("verification/Report.jl")
include("cli/main.jl")

end # module
