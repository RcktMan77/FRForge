# Publication / reproducibility snapshots.
#
# Freeze a method + reports + env lock into an immutable directory.
# Default verify is cheap (hashes + manifest); --rerun for full invent.

const SNAPSHOT_SCHEMA_VERSION = 1

"""
Known implementing source files for registered methods.
Convention: invent methods under `src/methods/<Name>.jl`.
Baselines may live under `src/capturing/`.
`--source` on freeze is an escape hatch only.
"""
const METHOD_SOURCE_MAP = Dict{String,Vector{String}}(
    "scaled_persson" => ["src/methods/ScaledPersson.jl"],
    "persson_av" => ["src/capturing/PerssonAV.jl", "src/capturing/PerssonAV2D.jl"],
    "null" => ["src/capturing/Interfaces.jl"],
)

"""
    resolve_method_sources(method_name; extra=String[]) -> Vector{String}

Registry-driven source discovery. Tries METHOD_SOURCE_MAP, then
`src/methods/<CamelOrName>.jl` heuristics. Appends `extra` paths.
"""
function resolve_method_sources(
    method_name::AbstractString;
    extra::AbstractVector{<:AbstractString}=String[],
)
    root = package_root()
    files = String[]
    if haskey(METHOD_SOURCE_MAP, String(method_name))
        append!(files, METHOD_SOURCE_MAP[String(method_name)])
    end
    # Heuristic: snake_case → ScaledPersson.jl style
    if isempty(files)
        parts = split(String(method_name), r"[_\-]+")
        camel = join(uppercasefirst.(parts))
        cand = joinpath("src", "methods", camel * ".jl")
        isfile(joinpath(root, cand)) && push!(files, cand)
        # also try exact name
        cand2 = joinpath("src", "methods", String(method_name) * ".jl")
        isfile(joinpath(root, cand2)) && push!(files, cand2)
    end
    for e in extra
        push!(files, String(e))
    end
    # unique preserve order
    seen = Set{String}()
    out = String[]
    for f in files
        f in seen && continue
        push!(seen, f)
        push!(out, f)
    end
    return out
end

function _file_sha256(path::AbstractString)
    data = read(path)
    return bytes2hex(SHA.sha256(data))
end

function _manifest_sha256()
    p = joinpath(package_root(), "Manifest.toml")
    isfile(p) || return "missing"
    return _file_sha256(p)
end

function _copy_file!(src::AbstractString, dest::AbstractString)
    mkpath(dirname(dest))
    cp(src, dest; force=true)
    return dest
end

"""
    freeze_snapshot(; method, baseline, method_report, baseline_report, compare,
                     out_root, git_ref, source_extra, log_path, append_log,
                     require_confirm, confirm_compare) -> String

Create an immutable snapshot directory. Returns absolute path.
Fails if zero source files resolve or if destination exists.

`require_confirm=true` hard-fails unless fine-mesh confirm evidence exists
(`frforge confirm` log entry or `confirm_compare` JSON). Default only warns.
"""
function freeze_snapshot(;
    method::AbstractString,
    baseline::AbstractString="persson_av",
    method_report::AbstractString,
    baseline_report::AbstractString,
    compare::Union{Nothing,AbstractString}=nothing,
    out_root::AbstractString=joinpath(package_root(), "results", "snapshots"),
    git_ref::AbstractString="",
    source_extra::AbstractVector{<:AbstractString}=String[],
    log_path::AbstractString=default_experiment_log_path(),
    append_log::Bool=true,
    date::Date=Dates.today(),
    require_confirm::Bool=false,
    confirm_compare::Union{Nothing,AbstractString}=nothing,
)
    root = package_root()
    isfile(method_report) || error("method_report not found: $method_report")
    isfile(baseline_report) || error("baseline_report not found: $baseline_report")

    conf_ok, conf_note = method_has_confirm_pass(
        method;
        log_path=log_path,
        confirm_compare=confirm_compare,
    )
    if !conf_ok
        msg = "Fine-mesh confirm not found for method=$method ($conf_note). " *
              "Paper-facing freezes should run `frforge confirm` first, then " *
              "`frforge snapshot freeze --require-confirm`."
        if require_confirm
            error(msg)
        else
            @warn msg
        end
    end

    sources = resolve_method_sources(method; extra=source_extra)
    # verify sources exist
    existing = String[]
    for s in sources
        p = isabspath(s) ? s : joinpath(root, s)
        if isfile(p)
            push!(existing, isabspath(s) ? relpath(p, root) : s)
        end
    end
    isempty(existing) &&
        error("No source files found for method \"$method\". Register in METHOD_SOURCE_MAP or pass --source.")

    sha = isempty(git_ref) ? git_commit_short() : String(git_ref)
    short = length(sha) > 8 ? sha[1:min(7, length(sha))] : sha
    stamp = Dates.format(date, dateformat"yyyymmdd")
    safe_method = replace(String(method), r"[^A-Za-z0-9_]+" => "_")
    snap_name = "$(safe_method)_$(stamp)_$(short)"
    snap_dir = joinpath(out_root, snap_name)
    isdir(snap_dir) && error("Snapshot already exists (immutable): $snap_dir")
    mkpath(snap_dir)
    mkpath(joinpath(snap_dir, "method"))
    mkpath(joinpath(snap_dir, "reports"))

    # copy sources
    rel_sources = String[]
    for s in existing
        src_path = joinpath(root, s)
        dest_name = basename(s)
        dest = joinpath(snap_dir, "method", dest_name)
        _copy_file!(src_path, dest)
        push!(rel_sources, joinpath("method", dest_name))
    end

    # copy reports
    reports = Dict{String,String}()
    mr_dest = joinpath(snap_dir, "reports", "method.json")
    _copy_file!(method_report, mr_dest)
    reports["method"] = "reports/method.json"
    br_dest = joinpath(snap_dir, "reports", "baseline.json")
    _copy_file!(baseline_report, br_dest)
    reports["baseline"] = "reports/baseline.json"
    if compare !== nothing && isfile(compare)
        cr_dest = joinpath(snap_dir, "reports", "compare.json")
        _copy_file!(compare, cr_dest)
        reports["compare"] = "reports/compare.json"
    end

    # extract experiment log entry (latest for method)
    entry_ids = String[]
    entry_md = ""
    if isfile(log_path)
        entries = parse_experiment_log(log_path)
        for e in reverse(entries)  # prefer latest
            if String(get(e, "method", "")) == String(method)
                push!(entry_ids, e["id"])
                entry_md = format_log_entry_text(e)
                break
            end
        end
    end
    open(joinpath(snap_dir, "experiment_entry.md"), "w") do io
        if isempty(entry_md)
            println(io, "# No matching experiment-log entry for method=$method")
        else
            write(io, entry_md)
        end
    end

    # load metrics from compare or method report
    primary = Dict{String,Any}()
    method_params = Dict{String,Any}()
    if haskey(reports, "compare")
        cmp = JSON.parsefile(joinpath(snap_dir, reports["compare"]))
        for k in (
            "composite",
            "candidate_status",
            "baseline_composite",
            "composite_margin",
            "tradeoff_ok",
        )
            haskey(cmp, k) && (primary[k] = cmp[k])
        end
        abs_s = get(cmp, "absolute_scores", Dict())
        if abs_s isa AbstractDict
            for k in ("order_preservation", "dissipation", "shock_quality", "robustness", "composite")
                haskey(abs_s, k) && (primary[k] = abs_s[k])
            end
        end
    end
    met_rep = JSON.parsefile(joinpath(snap_dir, reports["method"]))
    mp = get(met_rep, "method_params", Dict())
    mp isa AbstractDict && (method_params = Dict{String,Any}(String(k) => v for (k, v) in pairs(mp)))
    if !haskey(primary, "composite")
        scores = get(get(met_rep, "summary", Dict()), "scores", Dict())
        scores isa AbstractDict && (primary = merge(primary, Dict{String,Any}(String(k) => v for (k, v) in pairs(scores))))
    end

    # hashes of all packaged files
    hashes = Dict{String,String}()
    for (root_dir, _, files) in walkdir(snap_dir)
        for f in files
            p = joinpath(root_dir, f)
            rel = relpath(p, snap_dir)
            hashes[rel] = _file_sha256(p)
        end
    end
    # Store payload hashes only (hashes.json is the manifest of others).
    write_json_pretty(joinpath(snap_dir, "hashes.json"), hashes)

    snap = Dict{String,Any}(
        "schema_version" => SNAPSHOT_SCHEMA_VERSION,
        "method" => String(method),
        "baseline" => String(baseline),
        "git_ref" => sha,
        "created_at" => string(Dates.now()),
        "julia_version" => string(VERSION),
        "package_version" => package_version(),
        "manifest_sha256" => _manifest_sha256(),
        "scheme" => Dict{String,Any}(
            "points" => FROZEN_INVENT_SCHEME.points,
            "flux" => FROZEN_INVENT_SCHEME.flux,
            "time" => FROZEN_INVENT_SCHEME.time,
        ),
        "method_params" => method_params,
        "source_files" => rel_sources,
        "reports" => reports,
        "experiment_log_ids" => entry_ids,
        "primary_metrics" => primary,
        "re_run" => Dict{String,Any}(
            "command" => "frforge invent --method $method --baseline $baseline",
            "tol_rel" => 1e-10,
            "tol_abs" => 1e-12,
            "default_mode" => "cheap",
            "rerun_flag" => "--rerun",
        ),
        "snapshot_dir" => snap_name,
    )
    write_json_pretty(joinpath(snap_dir, "SNAPSHOT.json"), snap)

    # filled README from template
    readme = _fill_reproduce_readme(snap, snap_dir)
    open(joinpath(snap_dir, "README.md"), "w") do io
        write(io, readme)
    end

    if append_log
        _append_snapshot_log_entry!(
            log_path;
            method=method,
            git_ref=sha,
            snap_path=snap_dir,
            date=date,
        )
    end
    return abspath(snap_dir)
end

function _fill_reproduce_readme(snap::AbstractDict, snap_dir::AbstractString)
    io = IOBuffer()
    println(io, "# Reproducing method `$(snap["method"])`")
    println(io)
    println(io, "Snapshot: `$(basename(snap_dir))`  ")
    println(io, "Created: $(snap["created_at"])  ")
    println(io, "Git ref: `$(snap["git_ref"])`  ")
    println(io, "Julia: $(snap["julia_version"])  Package: $(snap["package_version"])")
    println(io)
    println(io, "## Environment")
    println(io)
    println(io, "1. Check out git ref `$(snap["git_ref"])` (or a commit that includes this method).")
    println(io, "2. `julia --project=. -e 'using Pkg; Pkg.instantiate()'`")
    println(io, "3. Manifest sha256 at freeze: `$(snap["manifest_sha256"])`")
    println(io)
    println(io, "## Cheap verify (default — no invent re-run)")
    println(io)
    println(io, "```bash")
    println(io, "frforge snapshot verify $(basename(snap_dir))")
    println(io, "# or with full path:")
    println(io, "frforge snapshot verify $snap_dir")
    println(io, "```")
    println(io)
    println(io, "## Full re-run (explicit)")
    println(io)
    println(io, "```bash")
    println(io, "frforge snapshot verify $snap_dir --rerun")
    println(io, "# equivalent invent:")
    println(io, snap["re_run"]["command"])
    println(io, "```")
    println(io)
    println(io, "## Tables from frozen JSON")
    println(io)
    println(io, "```bash")
    println(io, "frforge snapshot tables $snap_dir --out tables.md")
    println(io, "```")
    println(io)
    println(io, "## Primary metrics (frozen)")
    println(io)
    for (k, v) in sort(collect(pairs(snap["primary_metrics"])); by=x -> string(x[1]))
        println(io, "- **", k, ":** ", v)
    end
    println(io)
    println(io, "When to freeze: after invent short-list + fine-mesh confirm (+ robustness) — not after every invent run.")
    return String(take!(io))
end

"""
Short machine-friendly log append for snapshot_created (no full metrics dump).
"""
function _append_snapshot_log_entry!(
    log_path::AbstractString;
    method::AbstractString,
    git_ref::AbstractString,
    snap_path::AbstractString,
    date::Date=Dates.today(),
)
    id = make_entry_id(method; date=date, suffix="snapshot")
    # avoid id collision
    if isfile(log_path)
        existing = list_experiment_entry_ids(log_path)
        n = 1
        base = id
        while base in existing
            n += 1
            base = id * "-$n"
        end
        id = base
    end
    open(log_path, "a") do io
        println(io)
        println(io, "### ", id)
        println(io)
        println(io, "- **date:** ", Dates.format(date, dateformat"yyyy-mm-dd"))
        println(io, "- **method:** ", method)
        println(io, "- **baseline:** n/a")
        println(io, "- **hypothesis:** Snapshot freeze for reproducibility (not a new method invent).")
        println(io, "- **scheme:** points=", FROZEN_INVENT_SCHEME.points, ", flux=", FROZEN_INVENT_SCHEME.flux, ", time=", FROZEN_INVENT_SCHEME.time)
        println(io, "- **metrics:**")
        println(io, "  - candidate_status: snapshot_created")
        println(io, "- **lessons:** Use `frforge snapshot verify` (cheap) or `--rerun` for full invent.")
        println(io, "- **status:** snapshot_created")
        println(io, "- **artifacts:**")
        println(io, "  - snapshot: ", snap_path)
        println(io, "- **git_ref:** ", git_ref)
        println(io)
    end
    return id
end

"""
    verify_snapshot(snap_dir; rerun=false, tol_rel, tol_abs, require_git_ref=false) -> Dict

Default **cheap** path: hashes + SNAPSHOT.json + README presence.
With `rerun=true`, runs invent and compares metrics (slow — never required CI).
"""
function verify_snapshot(
    snap_dir::AbstractString;
    rerun::Bool=false,
    tol_rel::Real=1e-10,
    tol_abs::Real=1e-12,
    require_git_ref::Bool=false,
)
    snap_dir = abspath(snap_dir)
    isdir(snap_dir) || error("Not a directory: $snap_dir")
    snap_path = joinpath(snap_dir, "SNAPSHOT.json")
    isfile(snap_path) || error("Missing SNAPSHOT.json in $snap_dir")
    snap = JSON.parsefile(snap_path)
    hashes_path = joinpath(snap_dir, "hashes.json")
    isfile(hashes_path) || error("Missing hashes.json")
    hashes = JSON.parsefile(hashes_path)

    errors = String[]
    # hash check for files listed in hashes.json
    for (rel, expected) in pairs(hashes)
        p = joinpath(snap_dir, rel)
        if !isfile(p)
            push!(errors, "missing file $rel")
            continue
        end
        actual = _file_sha256(p)
        if actual != expected
            push!(errors, "hash mismatch $rel")
        end
    end
    isfile(joinpath(snap_dir, "README.md")) || push!(errors, "missing README.md")
    isfile(joinpath(snap_dir, "experiment_entry.md")) || push!(errors, "missing experiment_entry.md")

    if require_git_ref
        cur = git_commit_short()
        ref = string(get(snap, "git_ref", ""))
        if !isempty(ref) && ref != "unknown" && !startswith(cur, first(ref, min(7, length(ref)))) && cur != ref
            push!(errors, "git_ref mismatch: snapshot=$(ref) current=$(cur)")
        end
    end

    result = Dict{String,Any}(
        "ok" => isempty(errors),
        "mode" => rerun ? "rerun" : "cheap",
        "errors" => errors,
        "method" => get(snap, "method", ""),
        "snapshot_dir" => snap_dir,
        "primary_metrics" => get(snap, "primary_metrics", Dict()),
    )

    if rerun
        method = string(snap["method"])
        baseline = string(get(snap, "baseline", "persson_av"))
        # full invent — local/nightly only
        try
            _, _, cmp = invent_method(method; baseline=baseline, append_log=false)
            frozen = get(snap, "primary_metrics", Dict())
            abs_s = get(cmp, "absolute_scores", Dict())
            diffs = Dict{String,Any}()
            ok_metrics = true
            for k in ("composite", "order_preservation", "dissipation", "shock_quality")
                haskey(frozen, k) || continue
                fv = frozen[k]
                nv = get(abs_s, k, get(cmp, k, nothing))
                if fv isa Number && nv isa Number
                    if !isapprox(Float64(nv), Float64(fv); rtol=tol_rel, atol=tol_abs)
                        ok_metrics = false
                        diffs[k] = Dict("frozen" => fv, "rerun" => nv)
                    end
                end
            end
            result["rerun_ok"] = ok_metrics
            result["metric_diffs"] = diffs
            result["ok"] = result["ok"] && ok_metrics
            if !ok_metrics
                push!(errors, "metric mismatch on --rerun: $diffs")
                result["errors"] = errors
            end
        catch e
            push!(errors, "rerun failed: $e")
            result["ok"] = false
            result["errors"] = errors
        end
    end
    return result
end

"""
    snapshot_tables(snap_dir; out_md=nothing, out_csv=nothing) -> Dict

Regenerate comparison Markdown/CSV tables from frozen JSON reports.
"""
function snapshot_tables(
    snap_dir::AbstractString;
    out_md::Union{Nothing,AbstractString}=nothing,
    out_csv::Union{Nothing,AbstractString}=nothing,
)
    snap_dir = abspath(snap_dir)
    snap = JSON.parsefile(joinpath(snap_dir, "SNAPSHOT.json"))
    reports = get(snap, "reports", Dict())
    method_name = string(snap["method"])
    baseline_name = string(get(snap, "baseline", "persson_av"))

    rows = Dict{String,Any}[]
    if haskey(reports, "compare")
        cmp = JSON.parsefile(joinpath(snap_dir, reports["compare"]))
        abs_m = get(cmp, "absolute_scores", Dict())
        abs_b = get(cmp, "baseline_absolute_scores", Dict())
        for k in ("composite", "order_preservation", "dissipation", "shock_quality", "robustness")
            push!(
                rows,
                Dict{String,Any}(
                    "metric" => k,
                    "method" => get(abs_m, k, get(cmp, k, nothing)),
                    "baseline" => get(abs_b, k, nothing),
                ),
            )
        end
        push!(
            rows,
            Dict{String,Any}(
                "metric" => "candidate_status",
                "method" => get(cmp, "candidate_status", nothing),
                "baseline" => "baseline",
            ),
        )
    else
        met = get(snap, "primary_metrics", Dict())
        for (k, v) in pairs(met)
            push!(rows, Dict{String,Any}("metric" => k, "method" => v, "baseline" => nothing))
        end
    end

    md = IOBuffer()
    println(md, "# $(method_name) vs $(baseline_name)")
    println(md)
    println(md, "| metric | method | baseline |")
    println(md, "|--------|--------|----------|")
    for r in rows
        println(md, "| ", r["metric"], " | ", r["method"], " | ", r["baseline"], " |")
    end
    md_str = String(take!(md))

    csv = IOBuffer()
    println(csv, "metric,method,baseline")
    for r in rows
        println(csv, r["metric"], ",", r["method"], ",", r["baseline"])
    end
    csv_str = String(take!(csv))

    out_md !== nothing && (open(out_md, "w") do io
        write(io, md_str)
    end)
    out_csv !== nothing && (open(out_csv, "w") do io
        write(io, csv_str)
    end)

    return Dict{String,Any}(
        "markdown" => md_str,
        "csv" => csv_str,
        "rows" => rows,
        "method" => method_name,
        "baseline" => baseline_name,
    )
end
