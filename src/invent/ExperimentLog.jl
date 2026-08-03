# Persistent experiment log (laboratory notebook) for agent / human cumulative memory.
#
# Primary artifact: research/experiment_log.md (Markdown-primary)
# Optional index:  research/experiment_log.yaml
#
# Agents must read the log before proposing methods; invent appends after evaluation.

"""Frozen invent scheme for composite-score history (Phase 2+)."""
const FROZEN_INVENT_SCHEME = (
    points = "GL",
    flux = "Rusanov",
    time = "SSP-RK3",
)

const PROMISING_OR_HIGHER = Set([
    "promising",
    "accepted_candidate",
    "publication_grade",
])

const NARRATIVE_PLACEHOLDER =
    "[TODO: fill before shortlist/promotion — required for promising or higher]"

"""
    package_root() -> String

Repository root containing `Project.toml` and `research/`.
"""
function package_root()
    # src/invent/ExperimentLog.jl → src/invent → src → package root
    return dirname(dirname(@__DIR__))
end

"""
    default_experiment_log_path() -> String

Path to `research/experiment_log.md`.
"""
function default_experiment_log_path()
    return joinpath(package_root(), "research", "experiment_log.md")
end

"""
    default_experiment_log_yaml_path() -> String

Path to optional `research/experiment_log.yaml` index.
"""
function default_experiment_log_yaml_path()
    return joinpath(package_root(), "research", "experiment_log.yaml")
end

"""
    make_entry_id(method_name; date=today(), suffix="invent") -> String

Build a unique-ish entry id: `YYYYMMDD-method-suffix`.
"""
function make_entry_id(
    method_name::AbstractString;
    date::Date=Dates.today(),
    suffix::AbstractString="invent",
)
    safe = replace(String(method_name), r"[^A-Za-z0-9_]+" => "_")
    return string(Dates.format(date, dateformat"yyyymmdd"), "-", safe, "-", suffix)
end

"""
    narrative_required(candidate_status) -> Bool

True if hypothesis and lessons must be non-placeholder for this status.
"""
function narrative_required(candidate_status::AbstractString)
    return String(candidate_status) in PROMISING_OR_HIGHER
end

"""
    entry_from_invent(method_name, method_report, baseline_report, cmp;
                      hypothesis="", lessons="", strengths="", weaknesses="",
                      scheme=FROZEN_INVENT_SCHEME, git_ref="") -> Dict

Build a log entry dict from invent outputs.
"""
function entry_from_invent(
    method_name::AbstractString,
    method_report::AbstractDict,
    baseline_report::AbstractDict,
    cmp::AbstractDict;
    hypothesis::AbstractString="",
    lessons::AbstractString="",
    strengths::AbstractString="",
    weaknesses::AbstractString="",
    scheme=FROZEN_INVENT_SCHEME,
    git_ref::AbstractString="",
    artifacts::Union{Nothing,AbstractDict}=nothing,
    status::AbstractString="open",
    entry_id::Union{Nothing,AbstractString}=nothing,
    date::Date=Dates.today(),
)
    cand = String(get(cmp, "candidate_status", "unknown"))
    scores = get(cmp, "absolute_scores", get(get(method_report, "summary", Dict()), "scores", Dict()))
    hyp = String(hypothesis)
    les = String(lessons)
    if narrative_required(cand)
        isempty(strip(hyp)) && (hyp = NARRATIVE_PLACEHOLDER)
        isempty(strip(les)) && (les = NARRATIVE_PLACEHOLDER)
    else
        isempty(strip(hyp)) && (hyp = "Auto-appended invent run (edit to document hypothesis).")
        isempty(strip(les)) && (les = "Auto-appended invent run (edit to document lessons).")
    end

    arts = if artifacts !== nothing
        Dict{String,Any}(String(k) => v for (k, v) in pairs(artifacts))
    else
        Dict{String,Any}()
    end

    sch = if scheme isa NamedTuple
        Dict{String,Any}(
            "points" => string(scheme.points),
            "flux" => string(scheme.flux),
            "time" => string(scheme.time),
        )
    else
        Dict{String,Any}(
            "points" => string(get(scheme, "points", FROZEN_INVENT_SCHEME.points)),
            "flux" => string(get(scheme, "flux", FROZEN_INVENT_SCHEME.flux)),
            "time" => string(get(scheme, "time", FROZEN_INVENT_SCHEME.time)),
        )
    end

    return Dict{String,Any}(
        "id" => something(entry_id, make_entry_id(method_name; date=date)),
        "date" => string(date),
        "method" => String(method_name),
        "baseline" => string(something(get(cmp, "baseline_name", nothing), get(method_report, "baseline_name", "persson_av"))),
        "hypothesis" => hyp,
        "scheme" => sch,
        "metrics" => Dict{String,Any}(
            "candidate_status" => cand,
            "composite" => get(cmp, "composite", nothing),
            "baseline_composite" => get(cmp, "baseline_composite", nothing),
            "composite_margin" => get(cmp, "composite_margin", nothing),
            "order_preservation" => get(scores, "order_preservation", nothing),
            "dissipation" => get(scores, "dissipation", nothing),
            "shock_quality" => get(scores, "shock_quality", nothing),
            "robustness" => get(scores, "robustness", nothing),
            "tradeoff_ok" => get(cmp, "tradeoff_ok", nothing),
            "tradeoff_notes" => get(cmp, "tradeoff_notes", nothing),
        ),
        "strengths" => String(strengths),
        "weaknesses" => String(weaknesses),
        "lessons" => les,
        "status" => String(status),
        "artifacts" => arts,
        "git_ref" => String(git_ref),
        "narrative_complete" =>
            narrative_required(cand) ?
            (hyp != NARRATIVE_PLACEHOLDER && les != NARRATIVE_PLACEHOLDER) : true,
    )
end

"""
    format_entry_markdown(entry) -> String

Render one entry as a Markdown section (without trailing final blank line requirement).
"""
function format_entry_markdown(entry::AbstractDict)
    id = entry["id"]
    sch = entry["scheme"]
    metrics = entry["metrics"]
    io = IOBuffer()
    println(io, "### ", id)
    println(io)
    println(io, "- **date:** ", entry["date"])
    println(io, "- **method:** ", entry["method"])
    println(io, "- **baseline:** ", entry["baseline"])
    println(io, "- **hypothesis:** ", entry["hypothesis"])
    println(
        io,
        "- **scheme:** points=",
        get(sch, "points", "?"),
        ", flux=",
        get(sch, "flux", "?"),
        ", time=",
        get(sch, "time", "?"),
    )
    println(io, "- **metrics:**")
    for k in (
        "candidate_status",
        "composite",
        "baseline_composite",
        "composite_margin",
        "order_preservation",
        "dissipation",
        "shock_quality",
        "robustness",
        "tradeoff_ok",
        "tradeoff_notes",
    )
        if haskey(metrics, k) && metrics[k] !== nothing
            println(io, "  - ", k, ": ", metrics[k])
        end
    end
    str = get(entry, "strengths", "")
    isempty(str) || println(io, "- **strengths:** ", str)
    wkn = get(entry, "weaknesses", "")
    isempty(wkn) || println(io, "- **weaknesses:** ", wkn)
    println(io, "- **lessons:** ", get(entry, "lessons", ""))
    println(io, "- **status:** ", get(entry, "status", "open"))
    arts = get(entry, "artifacts", Dict())
    if arts isa AbstractDict && !isempty(arts)
        println(io, "- **artifacts:**")
        for (k, v) in sort(collect(pairs(arts)); by=x -> string(x[1]))
            println(io, "  - ", k, ": ", v)
        end
    else
        println(io, "- **artifacts:** none")
    end
    gr = get(entry, "git_ref", "")
    isempty(gr) || println(io, "- **git_ref:** ", gr)
    if get(entry, "narrative_complete", true) === false
        println(
            io,
            "- **note:** hypothesis/lessons still placeholders — required before shortlist for promising+",
        )
    end
    println(io)
    return String(take!(io))
end

"""
    append_experiment_entry!(path, entry; yaml_path=nothing) -> entry

Append a Markdown entry to the experiment log. Optionally append a one-line
YAML index note (simple append block) when `yaml_path` is provided and exists
or is creatable.

Returns the entry dict.
"""
function append_experiment_entry!(
    path::AbstractString,
    entry::AbstractDict;
    yaml_path::Union{Nothing,AbstractString}=nothing,
)
    mkpath(dirname(path))
    block = format_entry_markdown(entry)
    open(path, "a") do io
        # Ensure separation from previous content
        seekend(io)
        pos = position(io)
        if pos > 0
            # Ensure we start on a new section with blank line
            write(io, "\n")
        end
        write(io, block)
    end

    if yaml_path !== nothing
        _append_yaml_index_stub!(yaml_path, entry)
    end
    return entry
end

"""Minimal YAML index append (no full YAML rewrite dependency)."""
function _append_yaml_index_stub!(yaml_path::AbstractString, entry::AbstractDict)
    mkpath(dirname(yaml_path))
    metrics = get(entry, "metrics", Dict())
    open(yaml_path, "a") do io
        println(io)
        println(io, "  - id: ", entry["id"])
        println(io, "    date: \"", entry["date"], "\"")
        println(io, "    method: ", entry["method"])
        println(io, "    baseline: ", entry["baseline"])
        println(io, "    candidate_status: ", get(metrics, "candidate_status", "null"))
        println(io, "    status: ", get(entry, "status", "open"))
        if get(metrics, "composite", nothing) !== nothing
            println(io, "    composite: ", metrics["composite"])
        end
        arts = get(entry, "artifacts", Dict())
        if arts isa AbstractDict && !isempty(arts)
            println(io, "    artifacts:")
            for (k, v) in sort(collect(pairs(arts)); by=x -> string(x[1]))
                println(io, "      ", k, ": ", v)
            end
        end
    end
    return nothing
end

"""
    list_experiment_entry_ids(path) -> Vector{String}

Scan Markdown log for `### id` headings.
"""
function list_experiment_entry_ids(path::AbstractString)
    isfile(path) || return String[]
    ids = String[]
    for line in eachline(path)
        m = match(r"^###\s+(\S+)", line)
        if m !== nothing
            push!(ids, String(m.captures[1]))
        end
    end
    return ids
end

"""
    invent_append_log!(method_name, method_report, baseline_report, cmp;
                       log_path, yaml_path, artifacts, kwargs...) -> entry | nothing

Build entry from invent results and append to the experiment log.
"""
function invent_append_log!(
    method_name::AbstractString,
    method_report::AbstractDict,
    baseline_report::AbstractDict,
    cmp::AbstractDict;
    log_path::AbstractString=default_experiment_log_path(),
    yaml_path::Union{Nothing,AbstractString}=default_experiment_log_yaml_path(),
    artifacts::Union{Nothing,AbstractDict}=nothing,
    hypothesis::AbstractString="",
    lessons::AbstractString="",
    strengths::AbstractString="",
    weaknesses::AbstractString="",
    git_ref::AbstractString="",
    status::AbstractString="open",
)
    entry = entry_from_invent(
        method_name,
        method_report,
        baseline_report,
        cmp;
        hypothesis=hypothesis,
        lessons=lessons,
        strengths=strengths,
        weaknesses=weaknesses,
        artifacts=artifacts,
        git_ref=git_ref,
        status=status,
    )
    append_experiment_entry!(entry; path=log_path, yaml_path=yaml_path)
    return entry
end

# Convenience: keyword order matching common call style
function append_experiment_entry!(
    entry::AbstractDict;
    path::AbstractString=default_experiment_log_path(),
    yaml_path::Union{Nothing,AbstractString}=nothing,
)
    return append_experiment_entry!(path, entry; yaml_path=yaml_path)
end
