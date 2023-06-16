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
using LinearAlgebra
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

Convert value `z` to complex plane quadrant number.

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

Fill in `Vector` of candidate edges `E` that contain a phase change of 2 quadrants.

If `phasediffs` is not `nothing`, then push `|ΔQ|` of each edge to `phasediffs`.
This is useful for plotting.

Any root or pole is located at the point where the regions described by four different
quadrants meet. Since any triangulation of the four nodes located in the four different
quadrants requires at least one edge of ``|ΔQ| = 2``, then all such edges are potentially in
the vicinity of a root or pole.

`E` is not sorted.
"""
function candidateedges!(E, tess, pd=nothing)
    empty!(E)
    for edge in each_solid_edge(tess)
        _candidateedge!(E, tess, edge)
    end
end

function candidateedges!(E, tess, pd::PlotData)
    empty!(pd.phasediffs)
    for edge in each_solid_edge(tess)
        _candidateedge!(E, tess, edge)
        push!(pd.phasediffs, ΔQ)
    end
end

function _candidateedge!(E, tess, edge)
    a, b = get_point(tess, edge...)
    ΔQ = mod(getquadrant(a) - getquadrant(b), 4)  # phase difference
    if ΔQ == 2
        push!(E, edge)
    end
end

function selectedges!(selectE, tess, E, tolerance)
    empty!(selectE)
    for e in E
        d = distance(get_point(tess, e...))
        if d > tolerance
            push!(selectE, e)
        end
    end
end

"""
    addzone1node!(tess, p, q, tolerance)

Push the midpoint of `p` and `q` into `tess` if the distance between them is greater
than `tolerance`.
"""
function addzone1node!(tess, p, q, tolerance)
    if distance(p, q) > tolerance
        avgnode = (complex(p) + complex(q))/2
        add_point!(tess, QuadrantPoint(avgnode))
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
    splittriangles!(tess, unique_pts, tolerance, skinnytriangle)

Zone `1` triangles have more than one node in `edge_idxs`, whereas zone `2` triangles have
only a single node.
"""
function splittriangles!(tess, unique_pts, tolerance, skinnytriangle)
    triangles = Set{DT.triangle_type(tess)}()
    edges = Set{DT.edge_type(tess)}()
    for p1 in unique_pts
        adj2v = get_adjacent2vertex(tess, p1) 
        for (p2, p3) in adj2v
            sortedtri = DT.sort_triangle((p1, p2, p3))
            if sortedtri in triangles
                continue  # this triangle has already been visited
            end
            push!(triangles, sortedtri)
            p, q, r = get_point(tess, p1, p2, p3)

            if p2 in unique_pts || p3 in unique_pts
                # (p1, p2, p3) is a zone 1 triangle
                # Add a new node at the midpoint of each edge of (p1, p2, p3)
                if !(sort_edge(p1, p2) in edges)
                    addzone1node!(tess, p, q, tolerance)
                    push!(edges, sort_edge(p1, p2))
                end
                if !(sort_edge(p1, p3) in edges)
                    addzone1node!(tess, p, r, tolerance)
                    push!(edges, sort_edge(p1, p3))
                end
                if !(sort_edge(p2, p3) in edges)
                    addzone1node!(tess, q, r, tolerance)
                    push!(edges, sort_edge(p2, p3))
                end
            else
                # (p1, p2, p3) is a zone 2 triangle
                # Add a new node at the average of (p1, p2, p3) 
                addzone2node!(tess, p, q, r, skinnytriangle)
            end
        end
    end
end

"""
    contouredges(tess, E)

Find contour edges `C` from the edges of triangles containing at least one of candidate
edges `E`.
"""
function contouredges(tess, E)
    C = Set{DT.edge_type(tess)}()

    for e in E
        v = get_adjacent(tess, e)

        # If an edge occurs twice, that is because it is an edge shared by multiple triangles
        # and by definition is not an edge on the boundary of the candidate region.
        # Therefore, if an edge we come to is already in C, delete it.
        if (v, e[1]) in C
            delete!(C, (v, e[1]))
        elseif (e[1], v) in C
            # We need to check both edge orientations (a, b) and (b, a)
            delete!(C, (e[1], v))
        else
            push!(C, (v, e[1]))
        end

        if (v, e[2]) in C
            delete!(C, (v, e[2]))
        elseif (e[2], v) in C
            delete!(C, (e[2], v))
        else
            push!(C, (v, e[2]))
        end

        if e in C
            delete!(C, e)
        elseif reverse(e) in C
            delete!(C, reverse(e))
        else
            push!(C, e)
        end
    end

    return C
end

"""
    findnextpt(prevpt, refpt, nextedges)

Find the index of the next node in `nodes` as part of the candidate region boundary process.
The next one (after the reference) is picked from the fixed set of nodes.

Determined by the node that has the smallest difference in phase angle between the vector
from `prevnode` to `refnode` and the vector from the sample node to `refnode`. In other
words, it chooses the node that makes the least sharp among the possible paths. This traces
the outside closing contour around the candidate region.
"""
function findnextpt(tess, prevpt, refpt, nextedges)
    minphi = 2π + 1  # max diff of angles is 2π, so this is guaranteed larger
    c = first(nextedges)

    P = complex(get_point(tess, prevpt))
    S = complex(get_point(tess, refpt))
    for e in nextedges
        N = complex(get_point(tess, e[2]))
        SP = P - S
        SN = N - S

        phi = mod2pi(angle(SP) - angle(SN))

        if phi < minphi
            minphi = phi
            c = e
        end
    end

    return c
end

"""
    evaluateregions!(C, tess)

This function breaks C into the subsets Cᵏ for the closing contour surrounding the kth
candidate region.

Note: this function consumes `C`
"""
function evaluateregions!(C, tess)
    numregions = 1

    c = first(C)
    regions = [[c[1]]]
    refpt = c[2]
    delete!(C, c)

    while length(C) > 0
        nextedges = collect(Iterators.filter(v->v[1] == refpt, C))

        if isempty(nextedges)
            push!(regions[numregions], refpt)
            if length(C) > 0
                c = first(C)
                numregions += 1
                push!(regions, [c[1]])
                refpt = c[2]
                delete!(C, c) 
            end
        else
            if length(nextedges) > 1
                prevpt = regions[numregions][end]
                c = findnextpt(tess, prevpt, refpt, nextedges)
            else
                c = only(nextedges)
            end
            push!(regions[numregions], c[1])
            refpt = c[2]
            delete!(C, c)
        end
    end
    push!(regions[numregions], refpt)

    return regions
end


function evaluateregions(tess, C)
    regions = Vector{Vector{edge_type(tess)}}()
    assigned = Set{edge_type(tess)}()
    for c in C
        if !(c in assigned)
            for r in regions
                # XXX: `isdisjoint` isn't good enough - we need the countour edges ordered so we can apply DCAP
                if !isdisjoint(c, Iterators.flatten(r))
                    # `c` contains a vertex that is already in region `r`, so lets add it
                    # to the region
                    push!(r, c)
                end
            end
        end
    end
end

"""
    rootsandpoles(regions, quadrants)

Identify roots and poles of function based on `regions` and `quadrants`.
"""
function rootsandpoles(tess, regions)
    zroots = Vector{complex(number_type(tess))}()
    zpoles = similar(zroots)
    for r in regions
        quadrantsequence = [complex(get_point(tess, p)) for p in r]

        # Sign flip because `r` are in opposite order of Matlab?
        dq = -diff(quadrantsequence)
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
        z = sum(r)/length(r)

        if q > 0
            push!(zroots, convert(complexT, z))  # convert in case T isn't Float64
        elseif q < 0
            push!(zpoles, convert(complexT, z))
        end
    end

    return zroots, zpoles
end

"""
    tesselate!(tess, fcn, params, pd=nothing)

Label quadrants, identify candidate edges, and iteratively split triangles, returning
the tuple `(tess, E)`. Note that `tess` is modified in-place.
"""
function tesselate!(tess, fcn, params, pd::Union{Nothing,PlotData}=nothing)
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

        # Select candidate edges that are longer than the chosen tolerance
        selectedges!(selectE, tess, E, params.tolerance)
        isempty(selectE) && return tess, E

        # Get unique indices of points in `edges`
        unique_pts = Set(Iterators.flatten(selectE))

        # Refine tesselation
        # TODO: pass pd to store the point at every iteration
        splittriangles!(tess, unique_pts, params.tolerance, params.skinnytriangle)
    end

    # Assign quadrants to triangles split at end of last loop
    assignquadrants!(get_points(tess), fcn, params.multithreading)

    iteration >= params.maxiterations && @warn "$(iteration) iterations. `params.maxiterations` reached."
    num_vertices(tess) >= params.maxnodes && @warn "Tesselation has $(num_vertices(tess)) vertices. `params.maxnodes` reached."

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
    isempty(E) && return Vector{eltype(initial_mesh)}(), Vector{eltype(initial_mesh)}()

    C = contouredges(tess, E)
    regions = evaluateregions!(C)
    zroots, zpoles = rootsandpoles(regions, quadrants)

    return zroots, zpoles
end

end # module
