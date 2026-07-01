"""
    ComplexMesh{T<:Triangulation}

Struct that holds the mesh coordinates, Delaunay triangulation, function
values at the nodes of the mesh, and the random number generator for triangulation.

!!! warn

    DelaunayTriangulation.jl only supports triangulation of `Float64` meshes. 
"""
struct ComplexMesh{T<:Triangulation,R<:AbstractRNG}
    coords::Vector{NTuple{2,Float64}}
    tri::T
    fvals::Vector{ComplexF64}
    rng::R
    ComplexMesh{T,R}(c,t,f,r) where {T<:Triangulation,R<:AbstractRNG} =
        num_solid_vertices(t) == length(f) ?
        new(c,t,f,r) : error("Number of points is inconsistent")
end
ComplexMesh(coords, tri::T, rng::R) where {T<:Triangulation, R<:AbstractRNG} = ComplexMesh{T,R}(
    coords, tri, Vector{ComplexF64}(undef, length(coords)), rng
)

"""
    ComplexMesh(coords::Vector{NTuple{2,Float64}}; rng=Random.default_rng())

Create a `ComplexMesh` from `coords`. This form modifies `coords` in-place.

`rng` is passed to `DelaunayTriangulation.triangulate`.
"""
ComplexMesh(coords::Vector{NTuple{2,Float64}}; rng=Random.default_rng()) = ComplexMesh(coords, triangulate(coords; rng), rng)

"""
    ComplexMesh(coords::Vector{<:Number}; rng=Random.default_rng())

Attempt to promote `coords` to float type and split the real and imaginary parts.

This constructor builds a copy of `coords`, so it will not be modified in-place.
"""
ComplexMesh(coords::Vector{<:Number}; rng=Random.default_rng()) = ComplexMesh(reim.(float(coords)); rng)

function Base.show(io::IO, ::MIME"text/plain", m::ComplexMesh)
    println(io, "Complex Mesh – Delaunay Triangulation")
    println(io, "   Number of vertices: ", num_solid_vertices(m.tri))
    println(io, "   Number of triangles: ", num_solid_triangles(m.tri))
    println(io, "   Number of edges: ", num_solid_edges(m.tri))
    println(io, "   Has ghost triangles: ", DT.has_ghost_triangles(m.tri))
    println(io, "fvals: ", length(m.fvals), "-element ", typeof(m.fvals))
    println(io, "coords: ", length(m.coords), "-element ", typeof(m.coords))
    print(io, "rng: ", m.rng)
end

"Return the number of vertices in mesh `m`."
DT.num_solid_vertices(m::ComplexMesh) = DT.num_solid_vertices(m.tri)
Base.lastindex(m::ComplexMesh) = lastindex(m.coords)
Base.copy(m::ComplexMesh{T,R}) where {T,R} = ComplexMesh{T,R}(
    copy(m.coords), copy(m.tri), copy(m.fvals), m.rng
)

coord(m::ComplexMesh, i::Int) = m.coords[i]
coord(m::ComplexMesh, i::Vararg{Int, N}) where {N} = ntuple(j -> coord(m, i[j]), Val(N))
fval(m::ComplexMesh, i::Int) = m.fvals[i]
fval(m::ComplexMesh, i::Vararg{Int, N}) where {N} = ntuple(j -> fval(m, i[j]), Val(N))
fvals(m::ComplexMesh) = m.fvals
# using an fvals! in `evalfcn!` results in a runtime dispatch when `f` is an anonymous function
quadrant(m::ComplexMesh, i::Int) = quadrant(fval(m, i))
quadrant(m::ComplexMesh, i::Vararg{Int, N}) where {N} = ntuple(j -> quadrant(fval(m, i[j])), Val(N))

edge_length(m::ComplexMesh, e) = DT.edge_length(m.tri, e)::Float64
each_solid_edge(m::ComplexMesh) = DT.each_solid_edge(m.tri)
each_solid_triangle(m::ComplexMesh) = DT.each_solid_triangle(m.tri)
get_adjacent(m::ComplexMesh, e) = DT.get_adjacent(m.tri, e)
get_adjacent(m::ComplexMesh, a, b) = DT.get_adjacent(m.tri, a, b)
get_neighbours(m::ComplexMesh, M) = DT.get_neighbours(m.tri, M)
midpoint(m::ComplexMesh, u, v) = DT.midpoint(m.tri, u, v)

function addpoints!(m::ComplexMesh)
    n = num_solid_vertices(m.tri)
    nnew = length(m.coords)
    foreach(i -> add_point!(m.tri, i; m.rng), n+1:nnew)
    resize!(m.fvals, nnew)
    return nothing
end
