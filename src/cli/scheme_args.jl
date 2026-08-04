# Shared CLI scheme-axis argument helpers (behavior-preserving).
# Used by test/run/confirm parsers via scheme_from_cli_opts / help strings.

"""Add --points / --flux / --time to an ArgParseSettings table body via caller."""
const CLI_SCHEME_POINTS_HELP = "Solution points: gl (default) | gll"
const CLI_SCHEME_FLUX_HELP = "Numerical flux: rusanov (default) | hllc"
const CLI_SCHEME_TIME_HELP = "Time integrator: ssp_rk3 (default) | ssp_rk2"

"""Parse SchemeConfig from a CLI options dict with points/flux/time keys."""
function scheme_from_cli_opts(opts::AbstractDict)
    return parse_scheme(;
        points = get(opts, "points", "gl"),
        flux = get(opts, "flux", "rusanov"),
        time = get(opts, "time", "ssp_rk3"),
    )
end
