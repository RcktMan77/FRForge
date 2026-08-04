# Experiment-log analytics — read-only knowledge layer over research/experiment_log.md.
#
# Parses research/experiment_log.md into structured entries and produces
# summary / frontier / Pareto / lessons views for humans and agents.
# Never rewrites the authoritative Markdown log.

const FRONTIER_CANDIDATE = Set([
    "promising",
    "accepted_candidate",
    "publication_grade",
    "baseline",
    "confirmed",
])

const FRONTIER_LAB_STATUS = Set([
    "shortlisted",
    "publication_grade",
    "baseline",
    "confirmed",
])

"""Default composite margin for near-miss (pass_gates) inclusion on frontier."""
const NEAR_MISS_COMPOSITE_MARGIN = 0.02

# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

"""
    parse_experiment_log(path=default_experiment_log_path()) -> Vector{Dict{String,Any}}

Parse Markdown experiment log into normalized entry dicts.
Tolerant of missing optional fields, `~0.919` numerics, and platform entries.
"""
function parse_experiment_log(path::AbstractString = default_experiment_log_path())
    isfile(path) || return Dict{String,Any}[]
    text = read(path, String)
    return parse_experiment_log_text(text)
end

"""
    parse_experiment_log_text(text) -> Vector{Dict{String,Any}}

Parse from an in-memory Markdown string (also used for messy-fixture tests).
"""
function parse_experiment_log_text(text::AbstractString)
    entries = Dict{String,Any}[]
    # Split on ### headings that look like entry ids (not ## Schema etc.)
    parts = split(text, r"(?m)^###\s+")
    for (i, part) in enumerate(parts)
        i == 1 && continue  # preamble before first ###
        lines = split(part, '\n')
        isempty(lines) && continue
        id = strip(lines[1])
        # Skip non-entry ### headings (e.g. table docs — ids are YYYYMMDD-...)
        occursin(r"^\d{8}-", id) || continue
        body = join(lines[2:end], '\n')
        push!(entries, _parse_entry_body(id, body))
    end
    return entries
end

function _parse_entry_body(id::AbstractString, body::AbstractString)
    entry = Dict{String,Any}(
        "id" => String(id),
        "date" => "",
        "method" => "",
        "baseline" => "",
        "hypothesis" => "",
        "lessons" => "",
        "strengths" => "",
        "weaknesses" => "",
        "status" => "",
        "git_ref" => "",
        "scheme" => Dict{String,Any}("points" => nothing, "flux" => nothing, "time" => nothing),
        "metrics" => Dict{String,Any}(),
        "artifacts" => Dict{String,Any}(),
    )
    metrics = entry["metrics"]::Dict{String,Any}
    artifacts = entry["artifacts"]::Dict{String,Any}
    in_metrics = false
    in_artifacts = false

    for raw in split(body, '\n')
        line = rstrip(raw)
        isempty(strip(line)) && continue

        # Nested metrics: "  - key: value"
        if in_metrics
            m = match(r"^\s+-\s+([A-Za-z0-9_]+)\s*:\s*(.*)$", line)
            if m !== nothing
                key = String(m.captures[1])
                val = strip(m.captures[2])
                metrics[key] = _parse_metric_value(key, val)
                continue
            elseif startswith(strip(line), "- **")
                in_metrics = false
            else
                # prose continuation under metrics — ignore
                continue
            end
        end

        # Nested artifacts
        if in_artifacts
            m = match(r"^\s+-\s+([A-Za-z0-9_]+)\s*:\s*(.*)$", line)
            if m !== nothing
                artifacts[String(m.captures[1])] = strip(m.captures[2])
                continue
            elseif startswith(strip(line), "- **")
                in_artifacts = false
            else
                continue
            end
        end

        # Top-level bold field. Writer emits `- **field:** value` (colon inside bold).
        # Also accept `- **field**: value` and single-star variants.
        m = match(r"^-\s+\*\*([A-Za-z0-9_]+):?\*\*\s*:?\s*(.*)$", line)
        m === nothing && (m = match(r"^-\s+\*([A-Za-z0-9_]+):?\*\s*:?\s*(.*)$", line))
        m === nothing && continue
        field = lowercase(String(m.captures[1]))
        val = strip(m.captures[2])

        if field == "metrics"
            in_metrics = true
            in_artifacts = false
            if !isempty(val) && !startswith(val, "n/a")
                # inline metrics rare; store raw
                metrics["_inline"] = val
            end
        elseif field == "artifacts"
            in_artifacts = true
            in_metrics = false
            if occursin(r"(?i)none", val)
                # empty
            end
        elseif field == "scheme"
            entry["scheme"] = _parse_scheme_line(val)
        elseif field in ("date", "method", "baseline", "hypothesis", "lessons",
            "strengths", "weaknesses", "status", "git_ref", "note")
            entry[field == "note" ? "note" : field] = val
        else
            entry[field] = val
        end
    end

    # Normalize method for platform entries
    meth = String(get(entry, "method", ""))
    if occursin(r"(?i)platform|n/a|not a capturing", meth) || startswith(meth, "_(")
        entry["kind"] = "platform"
    else
        entry["kind"] = "method"
        # strip markdown italics leftovers
        entry["method"] = replace(meth, r"^_+|_+$" => "")
    end

    # candidate_status convenience at top level for filters
    if haskey(metrics, "candidate_status")
        entry["candidate_status"] = metrics["candidate_status"]
    end
    return entry
end

function _parse_scheme_line(val::AbstractString)
    sch = Dict{String,Any}("points" => nothing, "flux" => nothing, "time" => nothing)
    for (key, re) in (
        ("points", r"(?i)points\s*=\s*([A-Za-z0-9_\-]+)"),
        ("flux", r"(?i)flux\s*=\s*([A-Za-z0-9_\-]+)"),
        ("time", r"(?i)time\s*=\s*([A-Za-z0-9_\-]+)"),
    )
        m = match(re, val)
        m !== nothing && (sch[key] = String(m.captures[1]))
    end
    return sch
end

function _parse_metric_value(key::AbstractString, val::AbstractString)
    v = strip(val)
    isempty(v) && return nothing
    # boolean-ish
    if key in ("tradeoff_ok",)
        if occursin(r"(?i)^true", v)
            return true
        elseif occursin(r"(?i)^false", v)
            return false
        end
    end
    # strip leading ~ and take first float if present
    m = match(r"([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)", v)
    if m !== nothing && key in (
        "composite",
        "baseline_composite",
        "composite_margin",
        "order_preservation",
        "dissipation",
        "shock_quality",
        "robustness",
    )
        return tryparse(Float64, m.captures[1])
    end
    # candidate_status and free text: first token often enough
    if key == "candidate_status"
        tok = split(v)[1]
        return replace(tok, r"[^A-Za-z0-9_]" => "")
    end
    return v
end

"""
    get_experiment_entry(entries, id) -> Dict | nothing
"""
function get_experiment_entry(entries::AbstractVector, id::AbstractString)
    for e in entries
        String(get(e, "id", "")) == String(id) && return e
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Views
# ---------------------------------------------------------------------------

"""
    log_summary(entries) -> Dict

Counts by status / candidate_status, latest entry per method, open narrative TODOs.
"""
function log_summary(entries::AbstractVector)
    by_status = Dict{String,Int}()
    by_cand = Dict{String,Int}()
    latest = Dict{String,Dict{String,Any}}()
    todos = Dict{String,Any}[]
    confirm_failed = Dict{String,Any}[]

    for e in entries
        st = String(get(e, "status", "unknown"))
        by_status[st] = get(by_status, st, 0) + 1
        cs = string(
            get(
                e,
                "candidate_status",
                get(get(e, "metrics", Dict()), "candidate_status", "unknown"),
            ),
        )
        by_cand[cs] = get(by_cand, cs, 0) + 1

        meth = String(get(e, "method", ""))
        kind = String(get(e, "kind", "method"))
        if kind != "platform" && !isempty(meth) && meth != "n/a"
            # keep latest by date string then id
            if !haskey(latest, meth) ||
               string(e["date"], e["id"]) > string(latest[meth]["date"], latest[meth]["id"])
                latest[meth] = e
            end
        end

        hyp = String(get(e, "hypothesis", ""))
        les = String(get(e, "lessons", ""))
        if occursin("TODO", hyp) || occursin("TODO", les) || occursin(NARRATIVE_PLACEHOLDER, hyp) ||
           occursin(NARRATIVE_PLACEHOLDER, les)
            push!(
                todos,
                Dict{String,Any}(
                    "id" => e["id"],
                    "method" => meth,
                    "hypothesis_todo" =>
                        occursin("TODO", hyp) || occursin(NARRATIVE_PLACEHOLDER, hyp),
                    "lessons_todo" => occursin("TODO", les) || occursin(NARRATIVE_PLACEHOLDER, les),
                ),
            )
        end

        # Surface fine-mesh confirmation failures so agents do not retry blindly
        conf = string(get(get(e, "metrics", Dict()), "confirmation_status", ""))
        if st == "confirmation_failed" || conf == "confirmation_failed"
            push!(
                confirm_failed,
                Dict{String,Any}(
                    "id" => get(e, "id", ""),
                    "method" => meth,
                    "date" => get(e, "date", ""),
                    "caution" => "confirmation_failed — read lessons/artifacts before retrying fine-mesh confirm",
                    "mesh" => get(get(e, "metrics", Dict()), "mesh", nothing),
                ),
            )
        end
    end

    return Dict{String,Any}(
        "n_entries" => length(entries),
        "by_status" => by_status,
        "by_candidate_status" => by_cand,
        "latest_by_method" => Dict{String,Any}(
            k => Dict{String,Any}(
                "id" => v["id"],
                "date" => v["date"],
                "status" => get(v, "status", ""),
                "candidate_status" => get(v, "candidate_status", nothing),
                "composite" => get(get(v, "metrics", Dict()), "composite", nothing),
            ) for (k, v) in pairs(latest)
        ),
        "narrative_todos" => todos,
        "confirmation_failed" => confirm_failed,
    )
end

"""
    log_frontier(entries; near_margin=NEAR_MISS_COMPOSITE_MARGIN) -> Vector{Dict}

Core frontier (baseline + promising-class + shortlisted) plus near-miss pass_gates.
"""
function log_frontier(
    entries::AbstractVector;
    near_margin::Real = NEAR_MISS_COMPOSITE_MARGIN,
)
    # latest invent-like entry per method
    latest = Dict{String,Dict{String,Any}}()
    for e in entries
        kind = String(get(e, "kind", "method"))
        kind == "platform" && continue
        meth = String(get(e, "method", ""))
        (isempty(meth) || meth == "n/a") && continue
        if !haskey(latest, meth) ||
           string(e["date"], e["id"]) > string(latest[meth]["date"], latest[meth]["id"])
            latest[meth] = e
        end
    end

    # baseline composite reference
    base_comp = nothing
    for e in values(latest)
        cs = string(
            get(e, "candidate_status", get(get(e, "metrics", Dict()), "candidate_status", "")),
        )
        if cs == "baseline" || String(get(e, "status", "")) == "baseline"
            c = get(get(e, "metrics", Dict()), "composite", nothing)
            c isa Number && (base_comp = Float64(c))
        end
    end
    # best composite among invent methods for near-miss band
    best_comp = base_comp
    for e in values(latest)
        c = get(get(e, "metrics", Dict()), "composite", nothing)
        if c isa Number
            best_comp = best_comp === nothing ? Float64(c) : max(best_comp, Float64(c))
        end
    end

    rows = Dict{String,Any}[]
    for (meth, e) in pairs(latest)
        met = get(e, "metrics", Dict{String,Any}())
        cs = string(get(e, "candidate_status", get(met, "candidate_status", "")))
        st = String(get(e, "status", ""))
        comp = get(met, "composite", nothing)
        near = false
        core =
            cs in FRONTIER_CANDIDATE || st in FRONTIER_LAB_STATUS ||
            (cs == "baseline" || st == "baseline")
        if !core && cs == "pass_gates" && comp isa Number && best_comp !== nothing
            # competitive: within near_margin of best or baseline
            ref = best_comp
            if base_comp !== nothing
                ref = max(ref, base_comp)
            end
            if Float64(comp) >= ref - Float64(near_margin)
                near = true
            end
        end
        core || near || continue
        sch = get(e, "scheme", Dict())
        push!(
            rows,
            Dict{String,Any}(
                "method" => meth,
                "id" => e["id"],
                "date" => e["date"],
                "status" => st,
                "candidate_status" => cs,
                "near" => near && !core,
                "composite" => comp,
                "order_preservation" => get(met, "order_preservation", nothing),
                "dissipation" => get(met, "dissipation", nothing),
                "shock_quality" => get(met, "shock_quality", nothing),
                "robustness" => get(met, "robustness", nothing),
                "scheme" => sch,
            ),
        )
    end
    sort!(
        rows;
        by = r -> (
            -(r["near"] ? 0 : 1),
            -(r["composite"] isa Number ? Float64(r["composite"]) : -Inf),
            string(r["method"]),
        ),
    )
    return rows
end

"""
    log_pareto(entries) -> Vector{Dict}

Methods with numeric order / dissipation / shock scores; mark non-dominated.
"""
function log_pareto(entries::AbstractVector)
    latest = Dict{String,Dict{String,Any}}()
    for e in entries
        kind = String(get(e, "kind", "method"))
        kind == "platform" && continue
        meth = String(get(e, "method", ""))
        isempty(meth) && continue
        met = get(e, "metrics", Dict())
        o, d, s = get(met, "order_preservation", nothing),
        get(met, "dissipation", nothing),
        get(met, "shock_quality", nothing)
        (o isa Number && d isa Number && s isa Number) || continue
        if !haskey(latest, meth) ||
           string(e["date"], e["id"]) > string(latest[meth]["date"], latest[meth]["id"])
            latest[meth] = e
        end
    end
    rows = Dict{String,Any}[]
    for (meth, e) in pairs(latest)
        met = get(e, "metrics", Dict())
        push!(
            rows,
            Dict{String,Any}(
                "method" => meth,
                "id" => e["id"],
                "composite" => get(met, "composite", nothing),
                "order_preservation" => Float64(met["order_preservation"]),
                "dissipation" => Float64(met["dissipation"]),
                "shock_quality" => Float64(met["shock_quality"]),
                "pareto" => false,
            ),
        )
    end
    # non-dominated: maximize all three
    for i in eachindex(rows)
        ri = rows[i]
        dominated = false
        for j in eachindex(rows)
            i == j && continue
            rj = rows[j]
            if rj["order_preservation"] >= ri["order_preservation"] &&
               rj["dissipation"] >= ri["dissipation"] &&
               rj["shock_quality"] >= ri["shock_quality"] &&
               (
                   rj["order_preservation"] > ri["order_preservation"] ||
                   rj["dissipation"] > ri["dissipation"] ||
                   rj["shock_quality"] > ri["shock_quality"]
               )
                dominated = true
                break
            end
        end
        ri["pareto"] = !dominated
    end
    sort!(
        rows;
        by = r -> (-(r["composite"] isa Number ? Float64(r["composite"]) : -Inf), r["method"]),
    )
    return rows
end

"""
    log_lessons(entries; query=nothing) -> Vector{Dict}

Flatten lessons and weaknesses; optional case-insensitive substring filter.
"""
function log_lessons(entries::AbstractVector; query::Union{Nothing,AbstractString} = nothing)
    out = Dict{String,Any}[]
    q = query === nothing ? nothing : lowercase(String(query))
    for e in entries
        for (field, text) in
            (("lessons", get(e, "lessons", "")), ("weaknesses", get(e, "weaknesses", "")))
            t = strip(String(text))
            isempty(t) && continue
            if q !== nothing && !occursin(q, lowercase(t))
                continue
            end
            push!(
                out,
                Dict{String,Any}(
                    "id" => e["id"],
                    "method" => get(e, "method", ""),
                    "field" => field,
                    "text" => t,
                ),
            )
        end
    end
    return out
end

# ---------------------------------------------------------------------------
# Formatting (human text)
# ---------------------------------------------------------------------------

function format_log_summary_text(summary::AbstractDict)
    io = IOBuffer()
    println(io, "Experiment log summary  (n_entries=$(summary["n_entries"]))")
    println(io, "By status:")
    for (k, v) in sort(collect(pairs(summary["by_status"])); by = x -> x[1])
        println(io, "  ", k, ": ", v)
    end
    println(io, "By candidate_status:")
    for (k, v) in sort(collect(pairs(summary["by_candidate_status"])); by = x -> x[1])
        println(io, "  ", k, ": ", v)
    end
    println(io, "Latest per method:")
    for (m, info) in sort(collect(pairs(summary["latest_by_method"])); by = x -> x[1])
        println(
            io,
            "  ",
            m,
            "  id=",
            info["id"],
            "  status=",
            info["status"],
            "  cand=",
            info["candidate_status"],
            "  composite=",
            info["composite"],
        )
    end
    todos = summary["narrative_todos"]
    if !isempty(todos)
        println(io, "Narrative TODOs ($(length(todos))):")
        for t in todos
            println(io, "  ", t["id"], " method=", t["method"])
        end
    end
    fails = get(summary, "confirmation_failed", Any[])
    if !isempty(fails)
        println(io, "CAUTION — confirmation_failed ($(length(fails))):")
        for f in fails
            println(
                io,
                "  ",
                get(f, "id", "?"),
                " method=",
                get(f, "method", "?"),
                "  ",
                get(f, "caution", ""),
            )
        end
    end
    return String(take!(io))
end

function format_log_frontier_text(rows::AbstractVector)
    io = IOBuffer()
    println(io, "Frontier / near-miss methods")
    println(
        io,
        rpad("method", 18),
        rpad("cand", 18),
        rpad("near", 6),
        rpad("comp", 10),
        rpad("order", 8),
        rpad("dissip", 8),
        rpad("shock", 8),
        "id",
    )
    for r in rows
        println(
            io,
            rpad(string(r["method"]), 18),
            rpad(string(r["candidate_status"]), 18),
            rpad(r["near"] ? "yes" : "", 6),
            rpad(string(something(r["composite"], "")), 10),
            rpad(string(something(r["order_preservation"], "")), 8),
            rpad(string(something(r["dissipation"], "")), 8),
            rpad(string(something(r["shock_quality"], "")), 8),
            r["id"],
        )
    end
    return String(take!(io))
end

function format_log_pareto_text(rows::AbstractVector)
    io = IOBuffer()
    println(io, "Pareto-style order / dissipation / shock (pareto=non-dominated)")
    for r in rows
        println(
            io,
            rpad(string(r["method"]), 18),
            " order=",
            r["order_preservation"],
            " dissip=",
            r["dissipation"],
            " shock=",
            r["shock_quality"],
            " composite=",
            r["composite"],
            r["pareto"] ? "  [pareto]" : "",
        )
    end
    return String(take!(io))
end

function format_log_lessons_text(rows::AbstractVector)
    io = IOBuffer()
    println(io, "Lessons / weaknesses ($(length(rows)))")
    for r in rows
        println(io, "--- ", r["id"], " [", r["field"], "] method=", r["method"])
        println(io, r["text"])
        println(io)
    end
    return String(take!(io))
end

function format_log_entry_text(entry::AbstractDict)
    io = IOBuffer()
    println(io, "### ", entry["id"])
    for k in (
        "date",
        "method",
        "baseline",
        "status",
        "candidate_status",
        "hypothesis",
        "lessons",
        "strengths",
        "weaknesses",
        "git_ref",
    )
        haskey(entry, k) && !isempty(string(get(entry, k, ""))) &&
            println(io, k, ": ", entry[k])
    end
    sch = get(entry, "scheme", Dict())
    println(
        io,
        "scheme: points=",
        get(sch, "points", "?"),
        " flux=",
        get(sch, "flux", "?"),
        " time=",
        get(sch, "time", "?"),
    )
    met = get(entry, "metrics", Dict())
    if !isempty(met)
        println(io, "metrics:")
        for (k, v) in sort(collect(pairs(met)); by = x -> string(x[1]))
            println(io, "  ", k, ": ", v)
        end
    end
    arts = get(entry, "artifacts", Dict())
    if arts isa AbstractDict && !isempty(arts)
        println(io, "artifacts:")
        for (k, v) in sort(collect(pairs(arts)); by = x -> string(x[1]))
            println(io, "  ", k, ": ", v)
        end
    end
    return String(take!(io))
end
