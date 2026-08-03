# Normative key lists for JSON verification reports (schema_version = 1).
# See docs/design.md § JSON verification report schema.

const SCHEMA_VERSION = 1
const SCORING_FORMULA_VERSION = 1

const DEFAULT_SCORING_WEIGHTS = Dict{String,Float64}(
    "order_preservation" => 0.30,
    "dissipation" => 0.25,
    "shock_quality" => 0.25,
    "robustness" => 0.20,
)

"""Required top-level keys for any schema v1 report."""
const REQUIRED_TOP_LEVEL_KEYS = (
    "schema_version",
    "package",
    "package_version",
    "git_commit",
    "timestamp_utc",
    "julia_version",
    "command",
    "suite",
    "method_name",
    "method_params",
    "baseline_name",
    "overall_pass",
    "diverged",
    "nan_detected",
    "wall_time_sec",
    "scoring_weights",
    "scoring_formula_version",
    "hard_gate_failures",
    "cases",
    "summary",
)

"""Required keys under `summary`."""
const REQUIRED_SUMMARY_KEYS = (
    "n_cases",
    "n_passed",
    "n_failed",
    "scores",
)

"""Required keys under `summary.scores` when scores are present."""
const REQUIRED_SCORE_KEYS = (
    "order_preservation",
    "dissipation",
    "shock_quality",
    "robustness",
    "composite",
)

"""Common required fields on every case object."""
const REQUIRED_CASE_KEYS = (
    "name",
    "case_type",
    "equation",
    "p",
    "capturing_method",
    "pass",
    "diverged",
    "nan_detected",
    "conservation_residual",
    "conservation_pass",
    "conservation_metric",
    "positivity_ok",
    "wall_time_sec",
    "metrics",
)
