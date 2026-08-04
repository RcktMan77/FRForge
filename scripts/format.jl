#!/usr/bin/env julia
# Format src/ and test/ with the project's .JuliaFormatter.toml.
# Usage (from repo root):
#   julia --project=. -e 'using Pkg; Pkg.add("JuliaFormatter")'   # one-time / CI
#   julia --project=. scripts/format.jl
#   julia --project=. scripts/format.jl --check   # non-zero exit if changes needed

using JuliaFormatter

check_only = "--check" in ARGS
roots = ["src", "test"]
ok = true
for root in roots
    isdir(root) || continue
    # format returns true if already formatted
    already = format(root; verbose=false)
    if check_only && !already
        ok = false
        println(stderr, "Formatting needed under $root (run scripts/format.jl)")
    elseif !check_only
        println(already ? "OK already formatted: $root" : "Formatted: $root")
    end
end
check_only && !ok && exit(1)
check_only && println("Format check passed")
