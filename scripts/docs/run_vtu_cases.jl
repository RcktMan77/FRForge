#!/usr/bin/env julia
# Documentation-only VTU driver for publication figures.
# NOT used by invent or required CI.
#
# Usage (from repo root):
#   julia --project=. scripts/docs/run_vtu_cases.jl
#   julia --project=. scripts/docs/run_vtu_cases.jl --method persson_av
#   julia --project=. scripts/docs/run_vtu_cases.jl --method scaled_persson --outdir results/docs_vtu
#
# Writes high-order VTU with sensor + av diagnostics for:
#   1) 2D Riemann cfg6
#   2) Reduced Double-Mach-like

using FRForge
using ArgParse

function parse_cli(args)
    s = ArgParseSettings(
        description="Run documentation VTU cases (Riemann cfg6 + reduced DMR).",
        prog="run_vtu_cases.jl",
    )
    @add_arg_table! s begin
        "--method"
        help = "Capturing method (persson_av | scaled_persson recommended)"
        default = "persson_av"
        "--outdir"
        help = "Output directory for VTU files"
        default = "results/docs_vtu"
        "--tag"
        help = "Filename tag (default: method name)"
        default = ""
        "--riemann-n"
        arg_type = Int
        default = 48
        dest_name = "riemann_n"
        "--riemann-p"
        arg_type = Int
        default = 2
        dest_name = "riemann_p"
        "--riemann-t"
        arg_type = Float64
        default = 0.15
        dest_name = "riemann_t"
        "--dmr-nx"
        arg_type = Int
        default = 60
        dest_name = "dmr_nx"
        "--dmr-ny"
        arg_type = Int
        default = 20
        dest_name = "dmr_ny"
        "--dmr-p"
        arg_type = Int
        default = 1
        dest_name = "dmr_p"
        "--dmr-t"
        arg_type = Float64
        default = 0.08
        dest_name = "dmr_t"
    end
    return parse_args(args, s)
end

function main(args=ARGS)
    opts = parse_cli(args)
    method_name = opts["method"]
    method = get_capturing_method(method_name)
    tag = isempty(opts["tag"]) ? method_name : opts["tag"]
    outdir = opts["outdir"]
    mkpath(outdir)

    println("=== Documentation VTU driver ===")
    println("method=$method_name  outdir=$outdir  tag=$tag")
    println("scheme=GL+Rusanov+SSP-RK3 (frozen invent defaults; not reconfigured here)")

    # --- 1) Riemann cfg6 ---
    println("\n--- 2D Riemann cfg6 ---")
    c_r, state_r, eq_r = run_euler2d_riemann(;
        p=opts["riemann_p"],
        nx=opts["riemann_n"],
        ny=opts["riemann_n"],
        t_final=opts["riemann_t"],
        cfl=0.06,
        config=:cfg6,
        method=method,
        method_name=method_name,
    )
    path_r = joinpath(outdir, "riemann_cfg6_$(tag).vtu")
    if c_r["diverged"] || c_r["nan_detected"]
        @error "Riemann run failed" pass = c_r["pass"] diverged = c_r["diverged"]
    else
        write_vtu_high_order_with_capturing(path_r, state_r, eq_r, method)
        println("pass=$(c_r["pass"]) pos=$(c_r["positivity_ok"]) → $path_r")
    end

    # --- 2) Reduced Double Mach ---
    println("\n--- Double Mach (strength=:reduced) ---")
    c_d, state_d, eq_d = run_double_mach_reflection(;
        p=opts["dmr_p"],
        nx=opts["dmr_nx"],
        ny=opts["dmr_ny"],
        t_final=opts["dmr_t"],
        cfl=0.04,
        Lx=1.5,
        Ly=0.5,
        strength=:reduced,
        method=method,
        method_name=method_name,
        require_positivity=false,
    )
    path_d = joinpath(outdir, "double_mach_$(tag).vtu")
    if c_d["diverged"] || c_d["nan_detected"]
        @error "DMR run failed" pass = c_d["pass"] diverged = c_d["diverged"]
    else
        write_vtu_high_order_with_capturing(path_d, state_d, eq_d, method)
        println("pass=$(c_d["pass"]) pos=$(c_d["positivity_ok"]) → $path_d")
    end

    println("\nDone. Next: run ParaView scripts, e.g.")
    println("  pvpython scripts/docs/paraview/plot_2d_publication.py \\\\")
    println("    --vtu $path_r --outdir results/docs_figures --prefix riemann_cfg6_baseline")
    return 0
end

abspath(PROGRAM_FILE) == @__FILE__ && (main(); nothing)
