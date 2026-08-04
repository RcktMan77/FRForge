# Reusable residual workspaces (1D/2D) — avoid per-call zeros/similar without changing numerics.

"""Face-loop scratch for 2D residual (per-thread when residual threading is on)."""
struct FaceScratch2D{T}
    u_m::Vector{T}
    u_p::Vector{T}
    ug::Vector{T}
    fh::Vector{T}
end

FaceScratch2D(::Type{T}, Neq::Int) where {T} =
    FaceScratch2D{T}(
        Vector{T}(undef, Neq),
        Vector{T}(undef, Neq),
        Vector{T}(undef, Neq),
        Vector{T}(undef, Neq),
    )

"""2D residual workspace owned by SolutionState2D."""
mutable struct ResidualWorkspace2D{T}
    Np::Int
    Nel::Int
    Neq::Int
    u_work::Array{T,4}
    trW::Array{T,3}
    trE::Array{T,3}
    trS::Array{T,3}
    trN::Array{T,3}
    fhat_W::Array{T,3}
    fhat_E::Array{T,3}
    fhat_S::Array{T,3}
    fhat_N::Array{T,3}
    σ::Vector{T}
    face_pool::Vector{FaceScratch2D{T}}
    vol_pool::Vector{VolumeScratch2D{T}}
    # Persson sensor scratch (per-thread): U and V\\U temporaries
    sensor_U::Vector{Matrix{T}}
    sensor_Û::Vector{Matrix{T}}
    # BR0 AV pool (lazy size match)
    br0::Union{Nothing,NamedTuple}
end

function ResidualWorkspace2D(::Type{T}, Np::Int, Nel::Int, Neq::Int) where {T}
    nthr = max(Threads.nthreads(), 1)
    return ResidualWorkspace2D{T}(
        Np,
        Nel,
        Neq,
        zeros(T, Np, Np, Nel, Neq),
        zeros(T, Np, Neq, Nel),
        zeros(T, Np, Neq, Nel),
        zeros(T, Np, Neq, Nel),
        zeros(T, Np, Neq, Nel),
        zeros(T, Np, Neq, Nel),
        zeros(T, Np, Neq, Nel),
        zeros(T, Np, Neq, Nel),
        zeros(T, Np, Neq, Nel),
        zeros(T, Nel),
        [FaceScratch2D(T, Neq) for _ in 1:nthr],
        [VolumeScratch2D(T, Np, Neq) for _ in 1:nthr],
        [zeros(T, Np, Np) for _ in 1:nthr],
        [zeros(T, Np, Np) for _ in 1:nthr],
        nothing,
    )
end

@inline function sensor_scratch_for(ws::ResidualWorkspace2D)
    tid = min(Threads.threadid(), length(ws.sensor_U))
    return ws.sensor_U[tid], ws.sensor_Û[tid]
end

function ensure_residual_workspace!(state::SolutionState2D{T,Neq}) where {T,Neq}
    Np = n_points(state.ops)
    Nel = state.mesh.n_elements
    ws = state.residual_ws
    if ws === nothing || ws.Np != Np || ws.Nel != Nel || ws.Neq != Neq
        state.residual_ws = ResidualWorkspace2D(T, Np, Nel, Neq)
    end
    return state.residual_ws::ResidualWorkspace2D{T}
end

@inline function face_scratch_for(ws::ResidualWorkspace2D)
    pool = ws.face_pool
    length(pool) == 1 && return pool[1]
    return pool[min(Threads.threadid(), length(pool))]
end

"""1D residual workspace."""
mutable struct ResidualWorkspace1D{T}
    Np::Int
    Nel::Int
    Neq::Int
    u_work::Array{T,3}
    fL::Matrix{T}
    fR::Matrix{T}
    σ::Vector{T}
    traces::InterfaceTraces{T}
    u_m::Vector{T}
    u_p::Vector{T}
    f_vol::Matrix{T}       # (Np, Neq) physical flux scratch
    f_end_L::Vector{T}     # Neq discontinuous flux at ξ=-1
    f_end_R::Vector{T}     # Neq discontinuous flux at ξ=+1
    # 1D BR0 AV pool (lazy)
    br0::Union{Nothing,NamedTuple}
end

function ResidualWorkspace1D(::Type{T}, Np::Int, Nel::Int, Neq::Int) where {T}
    return ResidualWorkspace1D{T}(
        Np,
        Nel,
        Neq,
        zeros(T, Np, Nel, Neq),
        zeros(T, Nel, Neq),
        zeros(T, Nel, Neq),
        zeros(T, Nel),
        allocate_traces(Nel, Neq, T),
        Vector{T}(undef, Neq),
        Vector{T}(undef, Neq),
        zeros(T, Np, Neq),
        Vector{T}(undef, Neq),
        Vector{T}(undef, Neq),
        nothing,
    )
end

function ensure_residual_workspace!(state::SolutionState{T,Neq}) where {T,Neq}
    Np = n_points(state.ops)
    Nel = state.mesh.n_elements
    ws = state.residual_ws
    if ws === nothing || ws.Np != Np || ws.Nel != Nel || ws.Neq != Neq
        state.residual_ws = ResidualWorkspace1D(T, Np, Nel, Neq)
    end
    return state.residual_ws::ResidualWorkspace1D{T}
end

"""1D BR0 AV workspace; stored on ResidualWorkspace1D.br0."""
function ensure_br0_workspace_1d!(ws::ResidualWorkspace1D{T}, Np::Int, Nel::Int, Neq::Int) where {T}
    b = ws.br0
    if b !== nothing && b.Np == Np && b.Nel == Nel && b.Neq == Neq
        return b
    end
    b = (
        Np=Np,
        Nel=Nel,
        Neq=Neq,
        g=zeros(T, Np, Nel, Neq),
        uL=zeros(T, Nel, Neq),
        uR=zeros(T, Nel, Neq),
        gL=zeros(T, Nel, Neq),
        gR=zeros(T, Nel, Neq),
        ghat_L=zeros(T, Nel, Neq),
        ghat_R=zeros(T, Nel, Neq),
        ε=zeros(T, Nel),
    )
    ws.br0 = b
    return b
end

"""BR0 AV workspace (large arrays); stored on ResidualWorkspace2D.br0."""
function ensure_br0_workspace!(ws::ResidualWorkspace2D{T}, Np::Int, Nel::Int, Neq::Int) where {T}
    b = ws.br0
    if b !== nothing && b.Np == Np && b.Nel == Nel && b.Neq == Neq
        return b
    end
    b = (
        Np = Np,
        Nel = Nel,
        Neq = Neq,
        gx = zeros(T, Np, Np, Nel, Neq),
        gy = zeros(T, Np, Np, Nel, Neq),
        uW = zeros(T, Np, Nel, Neq),
        uE = zeros(T, Np, Nel, Neq),
        uS = zeros(T, Np, Nel, Neq),
        uN = zeros(T, Np, Nel, Neq),
        gxW = zeros(T, Np, Nel, Neq),
        gxE = zeros(T, Np, Nel, Neq),
        gyS = zeros(T, Np, Nel, Neq),
        gyN = zeros(T, Np, Nel, Neq),
        ghat_W = zeros(T, Np, Nel, Neq),
        ghat_E = zeros(T, Np, Nel, Neq),
        ghat_S = zeros(T, Np, Nel, Neq),
        ghat_N = zeros(T, Np, Nel, Neq),
    )
    ws.br0 = b
    return b
end
