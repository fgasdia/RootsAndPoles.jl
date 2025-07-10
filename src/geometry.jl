"""
    RPGeometry(tri::Triangulation)

This is a constructor for the [`RPGeometry`](@ref) struct, which holds the mesh and
associated data.

# Fields 
- `triangulation`: The underlying `Triangulation` from DelaunayTriangulation.jl.
- `triangulation_statistics`: The statistics of the triangulation. 
- `cv_volumes::Vector{Float64}`: A `Vector` of the volumes of each control volume.
- `triangle_props::Dict{NTuple{3,Int},TriangleProperties}`: A `Dict` mapping the indices of each triangle to its [`TriangleProperties`].
"""
struct RPGeometry{T}
    triangulation::T
    quadrants::Vector{Int}
    edges::Dict{}
end
function Base.show(io::IO, ::MIME"text/plain", geo::RPGeometry)
    nv = DelaunayTriangulation.num_solid_vertices(geo.triangulation_statistics)
    nt = DelaunayTriangulation.num_solid_triangles(geo.triangulation_statistics)
    print(io, "RPGeometry with $(nv) points and $(nt) triangles.")
end

get_quadrant(mesh::RPGeometry, i) = mesh.quadrants[i]
set_quadrant!(mesh::RPGeometry, i, q) = (mesh.quadrants[i] = q)

"""
    get_point(mesh, i)

Get the `i`th point in `mesh`.
"""
DT.get_point(mesh::RPGeometry, i) = DT.get_point(mesh.triangulation, i)

function RPGeometry(tri::Triangulation)
    nn = DelaunayTriangulation.num_points(tri)
    quadrants = zeros(nn)
    # triangle_props = Dict{NTuple{3,Int},TriangleProperties}()
    # for T in each_solid_triangle(tri)
    #     i, j, k = triangle_vertices(T)
    #     p, q, r = get_point(tri, i, j, k)
    #     triangle_props[triangle_vertices(T)] = TriangleProperties(shape_function_coefficients, ((m₁cx, m₁cy), (m₂cx, m₂cy), (m₃cx, m₃cy)), ((n₁x, n₁y), (n₂x, n₂y), (n₃x, n₃y)), (ℓ₁, ℓ₂, ℓ₃))
    # end
    return RPGeometry(tri, quadrants)
end
