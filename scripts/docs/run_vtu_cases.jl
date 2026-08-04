#!/usr/bin/env julia
# Documentation / presentation VTU driver (NOT invent, NOT required CI).
#
# Usage (from repo root):
#   # Presentation-quality meshes (default) — finer than CI-light
#   julia --project=. scripts/docs/run_vtu_cases.jl \
#     --method persson_av --outdir results/docs_vtu --tag baseline
#
#   # Explicit presets
#   julia --project=. scripts/docs/run_vtu_cases.jl --preset presentation --tag baseline
#   julia --project=. scripts/docs/run_vtu_cases.jl --preset quick        # coarser smoke
#
# Writes high-order VTU with sensor + av diagnostics for:
#   1) 2D Riemann cfg6
#   2) Reduced Double-Mach-like
#
# Frozen invent scheme GL+Rusanov+SSP-RK3 is not modified.

using FRForge
using ArgParse

"""Presentation vs quick (smoke) mesh presets. Polynomial degree and t_final match baseline docs."""
function apply_preset!(opts, preset::AbstractString)
    if preset == "presentation"
        # High-fidelity README/presentation meshes (much finer than CI-light).
        # Dense enough that HO tessellation + ResampleToImage hide element banding
        # on Cartesian Riemann shocks/contacts. Same p and t_final as prior docs baseline.
        # CI: Riemann ~16² p=1; DMR ~16×4.
        opts["riemann_n"] = 192
        opts["riemann_p"] = 2
        opts["riemann_t"] = 0.15
        opts["riemann_cfl"] = 0.035
        opts["dmr_nx"] = 280
        opts["dmr_ny"] = 100
        opts["dmr_p"] = 1
        opts["dmr_t"] = 0.08
        opts["dmr_cfl"] = 0.025
    elseif preset == "quick"
        opts["riemann_n"] = 32
        opts["riemann_p"] = 2
        opts["riemann_t"] = 0.12
        opts["riemann_cfl"] = 0.06
        opts["dmr_nx"] = 48
        opts["dmr_ny"] = 16
        opts["dmr_p"] = 1
        opts["dmr_t"] = 0.06
        opts["dmr_cfl"] = 0.04
    elseif preset == "ci"
        # Match CI-light gates (for smoke only)
        opts["riemann_n"] = 16
        opts["riemann_p"] = 1
        opts["riemann_t"] = 0.08
        opts["riemann_cfl"] = 0.08
        opts["dmr_nx"] = 16
        opts["dmr_ny"] = 4
        opts["dmr_p"] = 1
        opts["dmr_t"] = 0.03
        opts["dmr_cfl"] = 0.04
    else
        error("Unknown preset \"$preset\" (use presentation|quick|ci)")
    end
    return opts
end

function parse_cli(args)
    s = ArgParseSettings(
        description="Documentation VTU cases (Riemann cfg6 + reduced DMR).",
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
        "--preset"
        help = "Mesh preset: presentation (default) | quick | ci"
        default = "presentation"
        # Optional overrides (0 / NaN = use preset)
        "--riemann-n"
        arg_type = Int
        default = 0
        dest_name = "riemann_n"
        "--riemann-p"
        arg_type = Int
        default = -1
        dest_name = "riemann_p"
        "--riemann-t"
        arg_type = Float64
        default = NaN
        dest_name = "riemann_t"
        "--dmr-nx"
        arg_type = Int
        default = 0
        dest_name = "dmr_nx"
        "--dmr-ny"
        arg_type = Int
        default = 0
        dest_name = "dmr_ny"
        "--dmr-p"
        arg_type = Int
        default = -1
        dest_name = "dmr_p"
        "--dmr-t"
        arg_type = Float64
        default = NaN
        dest_name = "dmr_t"
        "--threads"
        help = "Opt-in residual threads (docs only; not for invent scores). Needs julia -t N."
        arg_type = Int
        default = 0
        "--progress-every"
        help = "Print progress every N SSP-RK3 steps (0 = auto when mesh is large)"
        arg_type = Int
        default = -1
        dest_name = "progress_every"
    end
    opts = parse_args(args, s)
    # materialize preset then apply explicit overrides
    base = Dict{String,Any}()
    apply_preset!(base, opts["preset"])
    for k in keys(base)
        haskey(opts, k) || (opts[k] = base[k])
    end
    # overrides only when user passed non-sentinel
    opts["riemann_n"] = opts["riemann_n"] > 0 ? opts["riemann_n"] : base["riemann_n"]
    opts["riemann_p"] = opts["riemann_p"] >= 0 ? opts["riemann_p"] : base["riemann_p"]
    opts["riemann_t"] = isnan(opts["riemann_t"]) ? base["riemann_t"] : opts["riemann_t"]
    opts["riemann_cfl"] = base["riemann_cfl"]
    opts["dmr_nx"] = opts["dmr_nx"] > 0 ? opts["dmr_nx"] : base["dmr_nx"]
    opts["dmr_ny"] = opts["dmr_ny"] > 0 ? opts["dmr_ny"] : base["dmr_ny"]
    opts["dmr_p"] = opts["dmr_p"] >= 0 ? opts["dmr_p"] : base["dmr_p"]
    opts["dmr_t"] = isnan(opts["dmr_t"]) ? base["dmr_t"] : opts["dmr_t"]
    opts["dmr_cfl"] = base["dmr_cfl"]
    return opts
end

function main(args=ARGS)
    opts = parse_cli(args)
    method_name = opts["method"]
    method = get_capturing_method(method_name)
    tag = isempty(opts["tag"]) ? method_name : opts["tag"]
    outdir = opts["outdir"]
    mkpath(outdir)

    # Docs-only threading: 0 → env FRFORGE_THREADS or 1; N → with_frforge_threads(N)
    thr = opts["threads"]
    thr = thr > 0 ? thr : frforge_thread_count()  # respect env if already set
    if opts["threads"] > 0
        thr = opts["threads"]
    end

    # Auto progress for large meshes (avoid "is it stuck?" anxiety)
    pe = opts["progress_every"]
    if pe < 0
        pe = (opts["riemann_n"] >= 48 || opts["dmr_nx"] >= 64) ? 200 : 0
    end

    println("=== Documentation VTU driver (preset=$(opts["preset"])) ===")
    println("method=$method_name  outdir=$outdir  tag=$tag")
    println("scheme=GL+Rusanov+SSP-RK3 (frozen invent defaults; not reconfigured)")
    println(
        "Riemann: p=$(opts["riemann_p"]) n=$(opts["riemann_n"])² t=$(opts["riemann_t"]) cfl=$(opts["riemann_cfl"])",
    )
    println(
        "DMR:     p=$(opts["dmr_p"]) nx=$(opts["dmr_nx"]) ny=$(opts["dmr_ny"]) t=$(opts["dmr_t"]) cfl=$(opts["dmr_cfl"])",
    )
    println(
        "Threaded residuals are for local documentation runs only; invent composite scores and promotion decisions always use the serial residual.",
    )

    return with_frforge_threads(max(1, thr)) do
        println(
            "residual threads: $(frforge_thread_count())  (julia nthreads=$(Threads.nthreads()); use julia -t N)",
        )
        pe > 0 && println("progress_every: $pe steps")

        # --- 1) Riemann cfg6 ---
        println("\n--- 2D Riemann cfg6 ---")
        t0 = time()
        c_r, state_r, eq_r = run_euler2d_riemann(;
            p=opts["riemann_p"],
            nx=opts["riemann_n"],
            ny=opts["riemann_n"],
            t_final=opts["riemann_t"],
            cfl=opts["riemann_cfl"],
            config=:cfg6,
            method=method,
            method_name=method_name,
            progress_every=pe,
            progress_label="riemann_cfg6",
        )
        path_r = joinpath(outdir, "riemann_cfg6_$(tag).vtu")
        if c_r["diverged"] || c_r["nan_detected"]
            @error "Riemann run failed" pass = c_r["pass"] diverged = c_r["diverged"]
        else
            write_vtu_high_order_with_capturing(path_r, state_r, eq_r, method)
            println(
                "pass=$(c_r["pass"]) pos=$(c_r["positivity_ok"]) wall=$(round(time()-t0;digits=1))s → $path_r",
            )
        end

        # --- 2) Reduced Double Mach ---
        println("\n--- Double Mach (strength=:reduced) ---")
        t0 = time()
        c_d, state_d, eq_d = run_double_mach_reflection(;
            p=opts["dmr_p"],
            nx=opts["dmr_nx"],
            ny=opts["dmr_ny"],
            t_final=opts["dmr_t"],
            cfl=opts["dmr_cfl"],
            Lx=1.5,
            Ly=0.5,
            strength=:reduced,
            method=method,
            method_name=method_name,
            require_positivity=false,
            progress_every=pe,
            progress_label="double_mach",
        )
        path_d = joinpath(outdir, "double_mach_$(tag).vtu")
        if c_d["diverged"] || c_d["nan_detected"]
            @error "DMR run failed" pass = c_d["pass"] diverged = c_d["diverged"]
        else
            write_vtu_high_order_with_capturing(path_d, state_d, eq_d, method)
            println(
                "pass=$(c_d["pass"]) pos=$(c_d["positivity_ok"]) wall=$(round(time()-t0;digits=1))s → $path_d",
            )
        end

        println("\nDone. Generate figures with:")
        println("  export PVPYTHON=/Applications/ParaView-6.1.0.app/Contents/bin/pvpython")
        println("  \$PVPYTHON scripts/docs/paraview/plot_2d_publication.py \\")
        println("    --vtu $path_r --outdir results/docs_figures --prefix riemann_cfg6_$(tag) --case riemann")
        println("  \$PVPYTHON scripts/docs/paraview/plot_2d_publication.py \\")
        println("    --vtu $path_d --outdir results/docs_figures --prefix double_mach_$(tag) --case dmr")
        return 0
    end
end

if abspath(PROGRAM_FILE) == abspath(@__FILE__)
    main()
end
