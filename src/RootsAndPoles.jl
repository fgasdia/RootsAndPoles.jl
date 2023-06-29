__precompile__(true)

"""
    RootsAndPoles.jl

`RootsAndPoles` is a Julia implementation of the Global complex Roots and Poles Finding
(GRPF) algorithm.

Matlab code is available under MIT license at https://github.com/PioKow/GRPF.

# References

[^1]: P. Kowalczyk, “Global complex roots and poles finding algorithm based on phase
    analysis for propagation and radiation problems,” IEEE Transactions on Antennas and
    Propagation, vol. 66, no. 12, pp. 7198–7205, Dec. 2018, doi: 10.1109/TAP.2018.2869213.
"""
module RootsAndPoles

import Base
import Base.Threads.@threads
using DelaunayTriangulation
const DT = DelaunayTriangulation

"""
    GRPFParams

Struct for holding values used by `RootsAndPoles` to stop iterating or split
Delaunay triangles.

Default values are provided by `GRPFParams()`.

# Fields

- `maxiterations::Int = 100`: the maximum number of refinement iterations before `grpf`
    returns.
- `maxnodes::Int = 500000`: the maximum number of Delaunay tessalation nodes before `grpf`
    returns.
- `skinnytriangle::Int = 3`: maximum ratio of the longest to shortest side length of
    Delaunay triangles before they are split during `grpf` refinement iterations.
- `tolerance::Float64 = 1e-9`: maximum allowed edge length of the tesselation defined in the
    `origcoords` domain before returning.
- `multithreading::Bool = false`: use `Threads.@threads` to run the user-provided function
    `fcn` across the `DelaunayTriangulation`.
"""
struct GRPFParams
    maxiterations::Int
    maxnodes::Int
    skinnytriangle::Int
    tolerance::Float64
    multithreading::Bool

    function GRPFParams(maxiterations, maxnodes, skinnytriangle, tolerance,
                        multithreading)
        new(maxiterations, maxnodes, skinnytriangle, tolerance, multithreading)
    end
end
GRPFParams() = GRPFParams(100, 500000, 3, 1e-9, false)

"""
    GRPFParams(tolerance, multithreading=false)

Convenience function for creating a `GRPFParams` object with the most important parameters,
`tolerance`, and `multithreading`.
"""
GRPFParams(tolerance, multithreading=false) =
    GRPFParams(100, 500000, 3, tolerance, multithreading)

function Base.isequal(a::GRPFParams, b::GRPFParams)
    for n in fieldnames(GRPFParams)
        isequal(getfield(a,n), getfield(b,n)) || return false
    end
    return true
end
Base.:(==)(a::GRPFParams, b::GRPFParams) = isequal(a,b)

struct PlotData
    phasediffs::Vector{Int}
end

# These files need the above structs defined
include("DelaunayTriangulationExtensions.jl")
include("utils.jl")
include("coordinate_domains.jl")

export rectangulardomain, diskdomain, grpf, PlotData, GRPFParams

"""
    quadrant(z)

Return complex plane quadrant of number `z`.

| Quadrant |       Phase       |
|:--------:|:-----------------:|
|    1     | 0 ≤ arg f < π/2   |
|    2     | π/2 ≤ arg f < π   |
|    3     | π ≤ arg f < 3π/2  |
|    4     | 3π/2 ≤ arg f < 2π |
"""
function quadrant(z)
    r, i = reim(z)
    if r > 0 && i >= 0
        return 1
    elseif r <= 0 && i > 0
        return 2
    elseif r < 0 && i <= 0
        return 3
    else
        # r >= 0 && i < 0
        return 4
    end
end

"""
    assignquadrants!(points, f, multithreading=false)

If the point quadrant is 0, evaluate function `f` for [`quadrant`](@ref) at each of `points`
and update each point in-place.
"""
function assignquadrants!(points, f, multithreading=false)
    if multithreading
        @threads for p in points
            if getquadrant(p) == 0
                Q = quadrant(f(complex(p)))
                setquadrant!(p, Q)
            end
        end
    else
        for p in points
            if getquadrant(p) == 0
                Q = quadrant(f(complex(p)))
                setquadrant!(p, Q)
            end
        end
    end
end

"""
    candidateedges!(E, tess, pd=nothing)

Empty candidate edges `E` and push edges from `tess` that contain a phase change of 2
quadrants.

If `phasediffs` is not `nothing`, then push `|ΔQ|` of each edge to `pd`. This is useful for
plotting.

Roots or poles are located where the regions described by four different quadrants meet.
Since any triangulation of the four nodes located in the four different quadrants requires
at least one edge of ``|ΔQ| = 2``, then all such edges are potentially in the vicinity of a
root or pole.

`E` is not sorted.
"""
function candidateedges!(E, tess, pd=nothing)
    empty!(E)
    for edge in each_solid_edge(tess)
        _candidateedge!(E, tess, edge)
    end
end

function candidateedges!(E, tess, pd::PlotData)
    empty!(E)
    empty!(pd.phasediffs)
    for edge in each_solid_edge(tess)
        _candidateedge!(E, tess, edge)
        push!(pd.phasediffs, ΔQ)
    end
end

function _candidateedge!(E, tess, edge)
    a, b = get_point(tess, edge[1], edge[2])
    ΔQ = mod(getquadrant(a) - getquadrant(b), 4)  # phase difference
    if ΔQ == 2
        push!(E, edge)
    end
end

function selectedges!(selectE, tess, E, tolerance)
    empty!(selectE)
    for e in E
        d = distance(get_point(tess, e[1], e[2]))
        if d > tolerance
            push!(selectE, e)
        end
    end
end

"""
    addzone1nodes!(tess, z1edges, tolerance)

Push the midpoint of each edge into `tess` if the distance between them is greater
than `tolerance` and the new vertex is unique.
"""
function addzone1nodes!(tess, z1edges, tolerance)
    newnodes = Vector{complex(DT.number_type(tess))}()
    for e in z1edges
        a, b = complex.(get_point(tess, e[1], e[2]))
        cn = (a + b)/2
        elength = distance(b, a)
        if elength > tolerance
            # Possible edge case where candidate `cn` is within 2eps from a new node?
            allseparate = true
            for n in newnodes
                dist = hypot(n - cn)
                if dist < 2*eps(DT.number_type(tess))
                    allseparate = false
                    break
                end
            end
            if allseparate
                push!(newnodes, cn)
                add_point!(tess, QuadrantPoint(cn))
            end
        end
    end
end

"""
    addzone2node!(tess, p, q, r, skinnytriangle)

Push the average of `p`, `q`, and `r` into `tess`.

`skinnytriangle` is the maximum allowed ratio of the longest to shortest side length.
"""
function addzone2node!(tess, p, q, r, skinnytriangle)
    l1 = distance(p, q)
    l2 = distance(p, r)
    l3 = distance(q, r)
    if max(l1,l2,l3)/min(l1,l2,l3) > skinnytriangle
        avgnode = (complex(p) + complex(q) + complex(r))/3
        add_point!(tess, QuadrantPoint(avgnode))
    end
end

"""
    splittriangles!(tess, unique_idxs, tolerance=GRPFParams().tolerance, skinnytriangle=GRPFParams().skinnytriangle)

Refine `tess` based on `unique_idxs` of select candidate edges.

See also [`addzone1node!`](@ref), [`addzone2node!`](@ref).
"""
function splittriangles!(tess, unique_idxs,
    tolerance=GRPFParams().tolerance, skinnytriangle=GRPFParams().skinnytriangle)

    # Determine which triangles are zone 1 and which are zone 2
    triangles = Dict{DT.triangle_type(tess),Int}()
    for w in unique_idxs
        adj2v = get_adjacent2vertex(tess, w)
        for (u, v) in adj2v
            if haskey(triangles, sorttriangle(u, v, w))
                triangles[sorttriangle(u, v, w)] += 1
            else
                triangles[sorttriangle(u, v, w)] = 1
            end
        end
    end

    z1edges = Vector{DT.edge_type(tess)}()
    for (t, c) in triangles
        u, v, w = indices(t)
        p, q, r = get_point(tess, u, v, w)
        if c == 1
            addzone2node!(tess, p, q, r, skinnytriangle)
        else
            # c > 1
            push!(z1edges, sortedge(u, v), sortedge(v, w), sortedge(w, u))
        end
    end
    sort!(z1edges)
    unique!(z1edges)
    addzone1nodes!(tess, z1edges, tolerance)
end

function triangleedges(tess, E)
    D = Vector{DT.edge_type(tess)}()
    for e in E
        # Get triangles sharing edge `e` - usually two for each edge
        v1 = get_neighbours(tess, e[1])
        v2 = get_neighbours(tess, e[2])
        vs = intersect(v1, v2)

        for v in vs
            tri = DT.construct_positively_oriented_triangle(tess, e[1], e[2], v)
            i, j, k = indices(tri)
            push!(D, (i, j), (j, k), (k, i))
        end
    end
    return unique(D)
end

"""
    contouredges(tess, E)

Return contour edges from the edges of triangles containing at least one of candidate
edges `E`.
"""
function contouredges(tess, E)
    C = Set{DT.edge_type(tess)}()

    for e in E 
        # Get triangles sharing edge `e`
        v1 = get_neighbours(tess, e[1])
        v2 = get_neighbours(tess, e[2])
        vs = intersect(v1, v2)
    
        # If an edge occurs twice, that is because it is an edge shared by multiple triangles
        # and by definition is not an edge on the boundary of the candidate region.
        # Therefore, if an edge we come to is already in C, delete it.
        for v in vs
            if (e[1], v) in C
                delete!(C, (e[1], v))
            elseif (v, e[1]) in C
                # We need to check both edge orientations (a, b) and (b, a)
                delete!(C, (v, e[1]))
            else
                push!(C, (e[1], v))
            end
    
            if (e[2], v) in C
                delete!(C, (e[2], v))
            elseif (v, e[2]) in C
                delete!(C, (v, e[2]))
            else
                push!(C, (e[2], v))
            end

            if e in C
                delete!(C, e)
            elseif (e[2], e[1]) in C
                delete!(C, (e[2], e[1]))
            else
                push!(C, e)
            end
        end
    end

    return C
end

"""
    findnextedge(tess, previdx, refidx, nextedges)

Return the edge in `nextedges` that traces out the candidate region boundary.

Determined by the node that has the smallest difference in phase angle between the vector
from `prevnode` to `refnode` and the vector from the sample node to `refnode`.
"""
function findnextedge(tess, previdx, refidx, nextedges)
    minphi = 2π + 1  # max diff of angles is 2π, so this is guaranteed larger
    c = first(nextedges)

    P = complex(get_point(tess, previdx))
    S = complex(get_point(tess, refidx))
    SP = P - S
    aSP = angle(SP)

    for e in nextedges
        N = complex(get_point(tess, e[2]))
        
        SN = N - S

        phi = mod2pi(aSP - angle(SN))

        if phi < minphi
            minphi = phi
            c = e
        end
    end

    return c
end

"""
    evaluateregions!(C, tess)

This function returns a `Vector` of `Vector`s of point indices that make up the contour of
the kth candidate region within `C`.

!!! note

    This function empties `C`.
"""
# TODO: we could rewrite this to work with C as a vector - we don't make use of it being a Set
function evaluateregions!(C, tess)
    c = first(C)
    regions = [[c[1]]]
    refidx = c[2]
    numregions = 1
    delete!(C, c)

    nextedges = Vector{DT.edge_type(tess)}()
    while length(C) > 0
        # Writing out the for loop avoids a Core.box closure issue with `refidx`
        # nextedges = [c for c in C if refidx in c]
        for c in C
            if c[1] == refidx
                push!(nextedges, c)
            end
        end

        if isempty(nextedges)
            push!(regions[numregions], refidx)

            c = first(C)
            numregions += 1
            push!(regions, [c[1]])
            refidx = c[2]
            delete!(C, c)
        else
            if length(nextedges) > 1
                previdx = last(regions[numregions])
                c = findnextedge(tess, previdx, refidx, nextedges)
            else
                c = only(nextedges)
            end
            push!(regions[numregions], c[1])
            refidx = c[2]
            delete!(C, c)
        end
        empty!(nextedges)
    end
    push!(regions[numregions], refidx)

    return regions
end


# function evaluateregions!(C, tess)
#     c = first(C)
#     regions = [[c[1]]]
#     refidx = c[2]
#     numregions = 1
#     delete!(C, c)

#     nextedges = Vector{DT.edge_type(tess)}()
#     while length(C) > 0
#         # Writing out the for loop avoids a Core.box closure issue with `refidx`
#         # nextedges = [c for c in C if refidx in c]
#         for c in C
#             if refidx in c
#                 push!(nextedges, c)
#             end
#         end

#         if isempty(nextedges)
#             # Contour has closed
#             push!(regions[numregions], refidx)
            
#             # Begin the next contour region
#             c = first(C)
#             numregions += 1
#             push!(regions, [c[1]])
#             refidx = c[2]
#             delete!(C, c)
#         else
#             if length(nextedges) > 1
#                 previdx = last(regions[numregions])
#                 c = findnextedge(tess, previdx, refidx, nextedges)
#                 ### TEMP XXX
#                 @info "length(nextedges) = $(length(nextedges))"
#             else
#                 c = only(nextedges)
#             end

#             if c[1] == refidx 
#                 push!(regions[numregions], c[1])
#                 refidx = c[2]
#             elseif c[2] == refidx
#                 push!(regions[numregions], c[2])
#                 refidx = c[1]
#             end
#             delete!(C, c)
#         end
#         empty!(nextedges)
#     end
#     push!(regions[numregions], refidx)

#     return regions
# end

"""
    rootsandpoles(tess, regions)

Return `Vector`s of roots and poles as identified from the collection `regions` of contour
point indices.

See also [`evaluateregions!`](@ref).
"""
function rootsandpoles(tess, regions)
    zroots = Vector{complex(DT.number_type(tess))}()
    zpoles = similar(zroots)
    for r in regions
        pts = get_point(tess, r...)
        quadrantsequence = [getquadrant(p) for p in pts]

        # Sign flip because `r` are in opposite order of Matlab?
        dq = diff(quadrantsequence)
        for i in eachindex(dq)
            if dq[i] == 3
                dq[i] = -1
            elseif dq[i] == -3
                dq[i] = 1
            elseif abs(dq[i]) == 2
                # ``|ΔQ| = 2`` is ambiguous; cannot tell whether phase increases or
                # decreases by two quadrants
                dq[i] = 0
            end
        end
        q = sum(dq)/4
        z = sum(complex, pts)/length(pts)

        if q > 0
            push!(zroots, z)  # convert in case T isn't Float64
        elseif q < 0
            push!(zpoles, z)
        end
    end

    return zroots, zpoles
end

"""
    tesselate!(tess, fcn, params=GRPFParams(), pd=nothing)

Label quadrants, identify candidate edges, and iteratively split triangles, returning
the tuple `(tess, E)` where `tess` is the refined tesselation and `E` is a collection of
edges that are candidates for being in the vicinity of a root or pole.
"""
function tesselate!(tess, fcn, params=GRPFParams(), pd::Union{Nothing,PlotData}=nothing)
    # It's more efficient to `empty!` the sets in the loop than it is to create and allocate
    # them from scratch
    E = Set{DT.edge_type(tess)}()  # edges
    selectE = Set{DT.edge_type(tess)}()

    iteration = 0
    while iteration < params.maxiterations && num_vertices(tess) < params.maxnodes
        iteration += 1

        # Evaluate the function `fcn` for the quadrant at each node
        assignquadrants!(get_points(tess), fcn, params.multithreading)

        # Determine candidate edges that may be near a root or pole
        # Candidate edges are those where the phase change |ΔQ| = 2
        candidateedges!(E, tess, pd)
        isempty(E) && return tess, E

        # Select candidate edges that are longer than `tolerance`
        selectedges!(selectE, tess, E, params.tolerance)
        isempty(selectE) && return tess, E

        # Get unique indices of points in `selectE`
        unique_idxs = Set(Iterators.flatten(selectE))

        # Refine tesselation
        # TODO: pass pd to store the point at every iteration
        splittriangles!(tess, unique_idxs, params.tolerance, params.skinnytriangle)
    end

    # Assign quadrants and candidate edges for triangles split at end of last loop
    assignquadrants!(get_points(tess), fcn, params.multithreading)
    # candidateedges!(E, tess, pd)

    iteration >= params.maxiterations && @info "$(iteration) iterations. `params.maxiterations` reached."
    num_vertices(tess) >= params.maxnodes && @info "Tesselation has $(num_vertices(tess)) vertices. `params.maxnodes` reached."

    return tess, E
end

"""
    grpf(fcn, initial_mesh, params=GRPFParams())

Return a vector `roots` and a vector `poles` of a single (complex) argument function
`fcn`.

Searches within a domain specified by the vector of complex `initial_mesh`.

# Examples
```jldoctest
julia> simplefcn(z) = (z - 1)*(z - im)^2*(z + 1)^3/(z + im)

julia> xb, xe = -2, 2

julia> yb, ye = -2, 2

julia> r = 0.1

julia> initial_mesh = rectangulardomain(complex(xb, yb), complex(xe, ye), r)

julia> roots, poles = grpf(simplefcn, initial_mesh);

julia> roots
3-element Array{Complex{Float64},1}:
    -0.9999999998508017 - 8.765385802810127e-11im
 1.3587683359414186e-10 + 1.0000000001862643im
     1.0000000002536367 + 1.0339757656912847e-26im

julia> poles
1-element Array{Complex{Float64},1}:
 -2.5363675604239815e-10 - 1.0000000002980232im
```
"""
function grpf(fcn, initial_mesh, params=GRPFParams())
    mesh_points = QuadrantPoints(QuadrantPoint.(initial_mesh))
    tess = triangulate(mesh_points)  # WARN: modifying `mesh_points` modifies `tess` in place
    
    tess, E = tesselate!(tess, fcn, params)

    # TODO: test for type stability when no roots or poles in domain of initial_mesh
    isempty(E) && return Vector{complex(DT.number_type(tess))}(), Vector{complex(DT.number_type(tess))}()

    C = contouredges(tess, E)
    regions = evaluateregions!(C, tess)
    zroots, zpoles = rootsandpoles(tess, regions)

    return zroots, zpoles
end

end # module
