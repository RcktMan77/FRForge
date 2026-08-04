# Agent invention surface: register methods under src/methods/ only.
#
# Layout (methods/):
#   Registry.jl        — list/describe/require + includes method files
#   ScaledPersson.jl   — example inventable method (composition demo)
#
# Workflow:
#   1. Create src/methods/MyMethod.jl implementing AbstractCapturingMethod hooks
#   2. include it from this file
#   3. register_method!("my_method", (; kwargs...) -> MyMethod(; kwargs...))
#   4. frforge invent --method my_method --baseline persson_av

"""List registered capturing method names (sorted)."""
function list_methods()
    return sort(collect(keys(METHOD_REGISTRY)))
end

"""
    describe_methods() -> String

Human-readable registry listing for CLI help.
"""
function describe_methods()
    io = IOBuffer()
    println(io, "Registered capturing methods:")
    for name in list_methods()
        println(io, "  - ", name)
    end
    return String(take!(io))
end

"""
    require_registered_method(name; role="method") -> String

Return `name` as `String`, or throw with the registered-method list.
Internal helper for CLI / invent / confirm entry points.
"""
function require_registered_method(name::AbstractString; role::AbstractString = "method")
    key = String(name)
    haskey(METHOD_REGISTRY, key) ||
        error("Unknown $role \"$key\". $(describe_methods())")
    return key
end

# Example inventable method (structural composition of Persson with different defaults)
include("ScaledPersson.jl")
register_method!("scaled_persson", (; kwargs...) -> ScaledPerssonMethod(; kwargs...))
