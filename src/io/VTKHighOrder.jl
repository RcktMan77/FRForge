# High-order VTK XML UnstructuredGrid writer (discontinuous Lagrange cells).
# ParaView ≥ 5.5 recommended for VTK_LAGRANGE_LINE / VTK_LAGRANGE_QUAD.

const VTK_LAGRANGE_LINE = 68
const VTK_LAGRANGE_QUAD = 70

"""
    vtk_lagrange_line_nodes(p; T=Float64) -> Vector{T}

Equispaced reference nodes on [-1,1] in **VTK Lagrange line order**:
1. vertices ξ=-1, ξ=+1
2. edge interiors left→right: ξ = -1 + 2k/p for k=1..p-1
"""
function vtk_lagrange_line_nodes(p::Int; T::Type=Float64)
    p >= 0 || throw(ArgumentError("p must be >= 0"))
    if p == 0
        return T[0]
    end
    nodes = Vector{T}(undef, p + 1)
    nodes[1] = -one(T)
    nodes[2] = one(T)
    for k in 1:(p - 1)
        nodes[2 + k] = -one(T) + T(2k) / T(p)
    end
    return nodes
end

"""
    gl_to_equi_interp(ops) -> Matrix

Interpolation matrix I such that u_equi = I * u_GL,
mapping GL solution-point values to VTK equispaced Lagrange nodes.
"""
function gl_to_equi_interp(ops::FROperators{T}) where {T}
    ξ_equi = vtk_lagrange_line_nodes(ops.p; T=T)
    return lagrange_basis_matrix(ops.ξ, ξ_equi)
end

"""Physical coordinates of VTK Lagrange nodes for element `e` (in VTK order)."""
function vtk_physical_nodes_1d(mesh::Mesh1D{T}, ops::FROperators{T}, e::Int) where {T}
    ξ = vtk_lagrange_line_nodes(ops.p; T=T)
    xL = mesh.x_vertices[e]
    xR = mesh.x_vertices[e + 1]
    mid = (xL + xR) / T(2)
    half = mesh.J[e]
    return mid .+ half .* ξ
end

"""
    write_vtu_high_order(path, state, eq; fields=:auto) -> path

Write discontinuous high-order VTU for a 1D `SolutionState`.

- Cell type: VTK_LAGRANGE_LINE (68)
- Topology: each element owns its own (p+1) nodes (face nodes duplicated)
- Fields:
  - scalar eq: `u` (conserved)
  - Euler: `rho`, `u`, `p` (primitives) and `rho`, `rho_u`, `E` (conserved)
"""
function write_vtu_high_order(
    path::AbstractString,
    state::SolutionState{T,Neq},
    eq;
    fields::Symbol=:auto,
) where {T,Neq}
    endswith(lowercase(path), ".vtu") ||
        @warn "VTK high-order writer produces .vtu; path does not end in .vtu" path

    mesh, ops = state.mesh, state.ops
    p = ops.p
    Np_vtk = p + 1
    Nel = mesh.n_elements
    n_points = Nel * Np_vtk
    n_cells = Nel

    Igl = gl_to_equi_interp(ops)

    # Build point coordinates (discontinuous) and interpolated fields
    coords = zeros(T, 3, n_points)  # x,y,z
    # Field storage: name => Vector of length n_points
    field_data = Dict{String,Vector{T}}()

    is_euler = eq isa Euler1D
    if fields === :auto
        if is_euler
            for name in ("rho", "u", "p", "rho_u", "E")
                field_data[name] = zeros(T, n_points)
            end
        else
            field_data["u"] = zeros(T, n_points)
        end
    end

    @inbounds for e in 1:Nel
        base = (e - 1) * Np_vtk  # 0-based start
        x_nodes = vtk_physical_nodes_1d(mesh, ops, e)
        for a in 1:Np_vtk
            pid = base + a
            coords[1, pid] = x_nodes[a]
            coords[2, pid] = zero(T)
            coords[3, pid] = zero(T)
        end

        # Interpolate each conserved component GL → equi
        U_equi = zeros(T, Np_vtk, Neq)
        for c in 1:Neq
            u_gl = @view state.u[:, e, c]
            U_equi[:, c] = Igl * u_gl
        end

        if is_euler
            for a in 1:Np_vtk
                pid = base + a
                ρ = U_equi[a, 1]
                ρu = U_equi[a, 2]
                E = U_equi[a, 3]
                vel = ρu / max(ρ, eps(T))
                # pressure with same EOS
                pres = pressure(eq, T[ρ, ρu, E])
                field_data["rho"][pid] = ρ
                field_data["u"][pid] = vel
                field_data["p"][pid] = pres
                field_data["rho_u"][pid] = ρu
                field_data["E"][pid] = E
            end
        else
            for a in 1:Np_vtk
                pid = base + a
                field_data["u"][pid] = U_equi[a, 1]
            end
        end
    end

    # Connectivity: sequential 0..p within each element block (VTK order matches storage)
    # Point layout per element is already VTK Lagrange order.
    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, "<?xml version=\"1.0\"?>")
        println(io, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">")
        println(io, "  <UnstructuredGrid>")
        println(io, "    <Piece NumberOfPoints=\"$n_points\" NumberOfCells=\"$n_cells\">")

        # Points
        println(io, "      <Points>")
        println(
            io,
            "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">",
        )
        for i in 1:n_points
            print(io, "          ", Float64(coords[1, i]), " ", Float64(coords[2, i]), " ", Float64(coords[3, i]))
            println(io)
        end
        println(io, "        </DataArray>")
        println(io, "      </Points>")

        # Cells
        println(io, "      <Cells>")
        println(io, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">")
        for e in 1:Nel
            base = (e - 1) * Np_vtk
            print(io, "          ")
            for a in 0:(Np_vtk - 1)
                print(io, base + a)
                a < Np_vtk - 1 && print(io, " ")
            end
            println(io)
        end
        println(io, "        </DataArray>")
        println(io, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">")
        print(io, "          ")
        for e in 1:Nel
            print(io, e * Np_vtk)
            e < Nel && print(io, " ")
        end
        println(io)
        println(io, "        </DataArray>")
        println(io, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">")
        print(io, "          ")
        for e in 1:Nel
            print(io, VTK_LAGRANGE_LINE)
            e < Nel && print(io, " ")
        end
        println(io)
        println(io, "        </DataArray>")
        println(io, "      </Cells>")

        # PointData
        scalar0 = is_euler ? "rho" : "u"
        println(io, "      <PointData Scalars=\"$scalar0\">")
        # Stable field order
        names = is_euler ? ("rho", "u", "p", "rho_u", "E") : ("u",)
        for name in names
            haskey(field_data, name) || continue
            println(io, "        <DataArray type=\"Float64\" Name=\"$name\" format=\"ascii\">")
            print(io, "          ")
            v = field_data[name]
            for i in 1:n_points
                print(io, Float64(v[i]))
                i < n_points && print(io, " ")
            end
            println(io)
            println(io, "        </DataArray>")
        end
        println(io, "      </PointData>")

        println(io, "    </Piece>")
        println(io, "  </UnstructuredGrid>")
        println(io, "</VTKFile>")
    end
    return path
end

"""
    vtk_point_counts_1d(n_elements, p) -> (n_points, n_cells, conn_len)

Oracle for discontinuous 1D Lagrange VTU sizing.
"""
function vtk_point_counts_1d(n_elements::Int, p::Int)
    return (n_elements * (p + 1), n_elements, p + 1)
end

"""
Parse a simple VTU we wrote: return NamedTuple with n_points, n_cells, types, connectivity sample.
"""
function parse_vtu_basic(path::AbstractString)
    txt = read(path, String)
    m_pts = match(r"NumberOfPoints=\"(\d+)\"", txt)
    m_cells = match(r"NumberOfCells=\"(\d+)\"", txt)
    n_points = m_pts === nothing ? -1 : parse(Int, m_pts.captures[1])
    n_cells = m_cells === nothing ? -1 : parse(Int, m_cells.captures[1])
    # types
    tm = match(r"Name=\"types\"[^>]*>(.*?)</DataArray>"s, txt)
    types = Int[]
    if tm !== nothing
        for tok in split(strip(tm.captures[1]))
            push!(types, parse(Int, tok))
        end
    end
    cm = match(r"Name=\"connectivity\"[^>]*>(.*?)</DataArray>"s, txt)
    conn = Int[]
    if cm !== nothing
        for tok in split(strip(cm.captures[1]))
            push!(conn, parse(Int, tok))
        end
    end
    return (n_points=n_points, n_cells=n_cells, types=types, connectivity=conn, text=txt)
end

# ---------------------------------------------------------------------------
# 2D Lagrange quadrilateral (VTK type 70)
# ---------------------------------------------------------------------------

"""
    vtk_lagrange_quad_nodes(p; T=Float64) -> Vector{NTuple{2,T}}

Reference (ξ,η) nodes on [-1,1]² in VTK Lagrange quad order:
corners, then edge interiors (bottom, right, top, left), then face interior.
"""
function vtk_lagrange_quad_nodes(p::Int; T::Type=Float64)
    p >= 1 || throw(ArgumentError("p >= 1 for quads"))
    nodes = NTuple{2,T}[]
    # corners: SW, SE, NE, NW
    push!(nodes, (-one(T), -one(T)))
    push!(nodes, (one(T), -one(T)))
    push!(nodes, (one(T), one(T)))
    push!(nodes, (-one(T), one(T)))
    if p >= 2
        # bottom edge interiors: ξ increasing
        for k in 1:(p - 1)
            ξ = -one(T) + T(2k) / T(p)
            push!(nodes, (ξ, -one(T)))
        end
        # right edge interiors: η increasing
        for k in 1:(p - 1)
            η = -one(T) + T(2k) / T(p)
            push!(nodes, (one(T), η))
        end
        # top edge interiors: ξ decreasing (from SE-ish to SW along top from vertex1 to vertex3? edge 2 is NE→NW)
        for k in 1:(p - 1)
            ξ = one(T) - T(2k) / T(p)
            push!(nodes, (ξ, one(T)))
        end
        # left edge interiors: η decreasing (NW→SW)
        for k in 1:(p - 1)
            η = one(T) - T(2k) / T(p)
            push!(nodes, (-one(T), η))
        end
        # face interior: (p-1)² tensor, ξ fastest, η outer, increasing
        for jj in 1:(p - 1), ii in 1:(p - 1)
            ξ = -one(T) + T(2ii) / T(p)
            η = -one(T) + T(2jj) / T(p)
            push!(nodes, (ξ, η))
        end
    end
    return nodes
end

function vtk_point_counts_2d(n_elements::Int, p::Int)
    n_per = (p + 1)^2
    return (n_elements * n_per, n_elements, n_per)
end

"""
    write_vtu_high_order(path, state::SolutionState2D, eq) -> path

Discontinuous high-order VTU with VTK_LAGRANGE_QUAD cells.
"""
function write_vtu_high_order(
    path::AbstractString,
    state::SolutionState2D{T,Neq},
    eq;
    fields::Symbol=:auto,
) where {T,Neq}
    mesh, ops = state.mesh, state.ops
    p = ops.p
    nodes_ref = vtk_lagrange_quad_nodes(p; T=T)
    n_per = length(nodes_ref)
    Nel = mesh.n_elements
    n_pts = Nel * n_per
    n_cells = Nel
    Np = length(ops.ξ)

    ξ_gl = ops.ξ
    ℓ_cache = [lagrange_at(ξ_gl, ξη[1]) for ξη in nodes_ref]
    m_cache = [lagrange_at(ξ_gl, ξη[2]) for ξη in nodes_ref]

    coords = zeros(T, 3, n_pts)
    is_euler = eq isa Euler2D
    field_data = Dict{String,Vector{T}}()
    if is_euler
        for name in ("rho", "u", "v", "p", "rho_u", "rho_v", "E")
            field_data[name] = zeros(T, n_pts)
        end
    else
        field_data["u"] = zeros(T, n_pts)
    end

    @inbounds for e in 1:Nel
        base = (e - 1) * n_per
        for a in 1:n_per
            ξ, η = nodes_ref[a]
            x, y = physical_xy(mesh, e, ξ, η)
            pid = base + a
            coords[1, pid] = x
            coords[2, pid] = y
            coords[3, pid] = zero(T)
            ℓξ = ℓ_cache[a]
            ℓη = m_cache[a]
            U = zeros(T, Neq)
            for c in 1:Neq
                s = zero(T)
                for jb in 1:Np, ib in 1:Np
                    s += state.u[ib, jb, e, c] * ℓξ[ib] * ℓη[jb]
                end
                U[c] = s
            end
            if is_euler
                ρ, u, v, pres = conserved_to_primitives(eq, U)
                field_data["rho"][pid] = ρ
                field_data["u"][pid] = u
                field_data["v"][pid] = v
                field_data["p"][pid] = pres
                field_data["rho_u"][pid] = U[2]
                field_data["rho_v"][pid] = U[3]
                field_data["E"][pid] = U[4]
            else
                field_data["u"][pid] = U[1]
            end
        end
    end

    mkpath(dirname(abspath(path)))
    open(path, "w") do io
        println(io, "<?xml version=\"1.0\"?>")
        println(io, "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">")
        println(io, "  <UnstructuredGrid>")
        println(io, "    <Piece NumberOfPoints=\"$n_pts\" NumberOfCells=\"$n_cells\">")
        println(io, "      <Points>")
        println(io, "        <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">")
        for i in 1:n_pts
            println(io, "          ", Float64(coords[1, i]), " ", Float64(coords[2, i]), " ", Float64(coords[3, i]))
        end
        println(io, "        </DataArray>")
        println(io, "      </Points>")
        println(io, "      <Cells>")
        println(io, "        <DataArray type=\"Int32\" Name=\"connectivity\" format=\"ascii\">")
        for e in 1:Nel
            base = (e - 1) * n_per
            print(io, "          ")
            for a in 0:(n_per - 1)
                print(io, base + a)
                a < n_per - 1 && print(io, " ")
            end
            println(io)
        end
        println(io, "        </DataArray>")
        println(io, "        <DataArray type=\"Int32\" Name=\"offsets\" format=\"ascii\">")
        print(io, "          ")
        for e in 1:Nel
            print(io, e * n_per)
            e < Nel && print(io, " ")
        end
        println(io)
        println(io, "        </DataArray>")
        println(io, "        <DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">")
        print(io, "          ")
        for e in 1:Nel
            print(io, VTK_LAGRANGE_QUAD)
            e < Nel && print(io, " ")
        end
        println(io)
        println(io, "        </DataArray>")
        println(io, "      </Cells>")
        scalar0 = is_euler ? "rho" : "u"
        println(io, "      <PointData Scalars=\"$scalar0\">")
        names = is_euler ? ("rho", "u", "v", "p", "rho_u", "rho_v", "E") : ("u",)
        for name in names
            println(io, "        <DataArray type=\"Float64\" Name=\"$name\" format=\"ascii\">")
            print(io, "          ")
            v = field_data[name]
            for i in 1:n_pts
                print(io, Float64(v[i]))
                i < n_pts && print(io, " ")
            end
            println(io)
            println(io, "        </DataArray>")
        end
        println(io, "      </PointData>")
        println(io, "    </Piece>")
        println(io, "  </UnstructuredGrid>")
        println(io, "</VTKFile>")
    end
    return path
end
