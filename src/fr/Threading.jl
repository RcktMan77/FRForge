# Opt-in residual threading for documentation / long local runs.
#
# Default is serial (bit-deterministic) for invent, scoring, confirm-for-promotion, CI.
# Enable only via FRFORGE_THREADS or with_frforge_threads / docs --threads.
#
# Threaded residual may differ by FP roundoff; never use for official invent composites.

# Uses LinearAlgebra + Base.Threads from parent module FRForge.

"""Task-local / scoped override for FRForge residual thread count (nothing → env default)."""
const _FRFORGE_THREADS_OVERRIDE = Ref{Union{Nothing,Int}}(nothing)

"""Previous BLAS thread count when we pin BLAS to 1 under multi-thread residual."""
const _BLAS_THREADS_SAVED = Ref{Union{Nothing,Int}}(nothing)

"""
    frforge_thread_count() -> Int

Effective residual worker count: `min(requested, Threads.nthreads())`, at least 1.
Request source (first match):
1. Scoped override from `with_frforge_threads`
2. Environment `FRFORGE_THREADS` (integer ≥ 1)
3. Default `1` (serial)

Note: `julia -t N` only makes workers available; FRForge still requires opt-in.
"""
function frforge_thread_count()
    ov = _FRFORGE_THREADS_OVERRIDE[]
    if ov !== nothing
        return max(1, min(Int(ov), nthreads()))
    end
    env = get(ENV, "FRFORGE_THREADS", "")
    if !isempty(env)
        n = tryparse(Int, env)
        if n !== nothing && n >= 1
            return max(1, min(n, nthreads()))
        end
    end
    return 1
end

"""True when residual kernels should use multi-threaded loops."""
residual_threading_enabled() = frforge_thread_count() > 1

"""
    with_frforge_threads(n::Integer) do ... end

Run `f` with residual thread request set to `n` (still capped by `Threads.nthreads()`).
When `n > 1`, pin OpenBLAS to 1 thread for the duration to avoid oversubscription.
Restores prior FRForge override and BLAS thread count on exit.
"""
function with_frforge_threads(f::Function, n::Integer)
    n = max(1, Int(n))
    prev = _FRFORGE_THREADS_OVERRIDE[]
    _FRFORGE_THREADS_OVERRIDE[] = n
    blas_pinned = false
    if n > 1 && nthreads() > 1
        try
            _BLAS_THREADS_SAVED[] = BLAS.get_num_threads()
            BLAS.set_num_threads(1)
            blas_pinned = true
        catch
            # BLAS backend may not support query/set
        end
    end
    try
        return f()
    finally
        _FRFORGE_THREADS_OVERRIDE[] = prev
        if blas_pinned
            saved = _BLAS_THREADS_SAVED[]
            _BLAS_THREADS_SAVED[] = nothing
            if saved !== nothing
                try
                    BLAS.set_num_threads(saved)
                catch
                end
            end
        end
    end
end

"""Force serial residual for invent / score / official confirm."""
with_serial_residual(f::Function) = with_frforge_threads(f, 1)

"""
    foreach_element(f, Nel)

Call `f(e)` for each element `e ∈ 1:Nel`. Serial when `frforge_thread_count()==1`
(preserves loop order for bit-determinism). Otherwise static chunks over up to
`frforge_thread_count()` workers.
"""
function foreach_element(f::F, Nel::Int) where {F}
    nt = frforge_thread_count()
    if nt <= 1 || nthreads() == 1 || Nel < 2
        @inbounds for e in 1:Nel
            f(e)
        end
        return nothing
    end
    nt = min(nt, nthreads())
    chunk = cld(Nel, nt)
    @sync for tid in 1:nt
        Threads.@spawn begin
            e0 = (tid - 1) * chunk + 1
            e1 = min(tid * chunk, Nel)
            @inbounds for e in e0:e1
                f(e)
            end
        end
    end
    return nothing
end

"""Per-thread volume residual scratch (allocated once per residual call)."""
struct VolumeScratch2D{T}
    Fx::Vector{T}
    Gy::Vector{T}
    Ft::Array{T,3}
    Gt::Array{T,3}
end

function VolumeScratch2D(::Type{T}, Np::Int, Neq::Int) where {T}
    return VolumeScratch2D{T}(
        Vector{T}(undef, Neq),
        Vector{T}(undef, Neq),
        zeros(T, Np, Np, Neq),
        zeros(T, Np, Np, Neq),
    )
end

"""
    make_volume_scratch_pool(T, Np, Neq) -> Vector{VolumeScratch2D}

Pool length = `max(frforge_thread_count(), 1)` (or `nthreads()` when threaded).
Index with `Threads.threadid()` when multi-threaded; use `pool[1]` when serial.
"""
function make_volume_scratch_pool(::Type{T}, Np::Int, Neq::Int) where {T}
    nt = residual_threading_enabled() ? max(frforge_thread_count(), nthreads()) : 1
    return [VolumeScratch2D(T, Np, Neq) for _ in 1:nt]
end

@inline function volume_scratch_for(pool::Vector{<:VolumeScratch2D})
    if length(pool) == 1
        return pool[1]
    end
    tid = Threads.threadid()
    return pool[min(tid, length(pool))]
end

"""Per-thread 1D work vector for local DD AV."""
function make_vector_scratch_pool(::Type{T}, n::Int) where {T}
    nt = residual_threading_enabled() ? max(frforge_thread_count(), nthreads()) : 1
    return [Vector{T}(undef, n) for _ in 1:nt]
end

@inline function vector_scratch_for(pool::Vector{<:AbstractVector})
    if length(pool) == 1
        return pool[1]
    end
    tid = Threads.threadid()
    return pool[min(tid, length(pool))]
end
