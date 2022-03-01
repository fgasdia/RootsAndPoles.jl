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

#==
NOTE: Some variable conversions from the original GRPF papers to this code:

| Paper | Code |
|-------|------|
|   𝓔   |   E  |
|   𝐶   |   C  |
|   ϕ   |  phi |
==#

import Base
import Base.Threads.@threads
using LinearAlgebra
using VoronoiDelaunay
import VoronoiDelaunay: getx, gety, DelaunayEdge, DelaunayTriangle

# NOTE: `max_coord` and `min_coord` are provided by `VoronoiDelaunay.jl`
# We are even more conservative than going from `max_coord` to `min_coords` because it is
# possible to run into floating point issues at the very limits
const MAXCOORD = nextfloat(max_coord, -10)
const MINCOORD = nextfloat(min_coord, 10)

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
- `tess_sizehint::Int = 5000`: provide a size hint to the total number of expected nodes in
    the Delaunay tesselation. Setting this number approximately correct can improve
    performance.
- `tolerance::Float64 = 1e-9`: maximum allowed edge length of the tesselation defined in the
    `origcoords` domain before returning.
- `multithreading::Bool = false`: use `Threads.@threads` to run the user-provided function
    `fcn` across the `DelaunayTriangulation`.
"""
struct GRPFParams
    maxiterations::Int
    maxnodes::Int
    skinnytriangle::Int
    tess_sizehint::Int
    tolerance::Float64
    multithreading::Bool

    function GRPFParams(maxiterations, maxnodes, skinnytriangle, tess_sizehint, tolerance,
                        multithreading)
        tess_sizehint > maxnodes && @warn "GRPFParams `tess_sizehint` is greater than `maxnodes`"
        new(maxiterations, maxnodes, skinnytriangle, tess_sizehint, tolerance, multithreading)
    end
end
GRPFParams() = GRPFParams(100, 500000, 3, 5000, 1e-9, false)

"""
    GRPFParams(tess_sizehint, tolerance, multithreading=false)

Convenience function for creating a `GRPFParams` object with the most important parameters,
`tess_sizehint`, `tolerance`, and `multithreading`.
"""
GRPFParams(tess_sizehint, tolerance, multithreading=false) =
    GRPFParams(100, 500000, 3, tess_sizehint, tolerance, multithreading)

function Base.isequal(a::GRPFParams, b::GRPFParams)
    for n in fieldnames(GRPFParams)
        isequal(getfield(a,n), getfield(b,n)) || return false
    end
    return true
end
Base.:(==)(a::GRPFParams, b::GRPFParams) = isequal(a,b)

struct PlotData end

"""
    Geometry2Function{T}

Store conversion coefficients from the `VoronoiDelaunay` domain to the `origcoords`
function domain.

`Geometry2Function` instances can be called as a function, e.g.
`z_functiondomain = g2f(z_voronoidelaunay)`.

`Geometry2Function` structs are created internally to `RootsAndPoles` with the
appropriate scaling parameters but are returned when `grpf` is called with the `PlotData()`
argument.
"""
struct Geometry2Function{T}
    rmin::T
    rmax::T
    imin::T
    imax::T
end
(f::Geometry2Function)(z) = geom2fcn(z, f.rmin, f.rmax, f.imin, f.imax)
(f::Geometry2Function)(x, y) = geom2fcn(x, y, f.rmin, f.rmax, f.imin, f.imax)
Base.eltype(f::Geometry2Function{T}) where T = T

"""
    ScaledFunction{F,G<:Geometry2Function}

Store conversion coefficients from the `VoronoiDelaunay` domain to the `origcoords` domain
so that a `ScaledFunction` can be called in place of the original function when providing an
argument in the `VoronoiDelaunay` domain.
"""
struct ScaledFunction{F,G<:Geometry2Function}
    fcn::F
    g2f::G
end
(f::ScaledFunction)(z) = f.fcn(f.g2f(z))

# These files need the above structs defined
include("VoronoiDelaunayExtensions.jl")
include("utils.jl")
include("coordinate_domains.jl")

export rectangulardomain, diskdomain, grpf, PlotData, getplotdata, GRPFParams

"""
    quadrant(val)

Convert complex function value `val` to quadrant number.

| Quadrant |       Phase       |
|:--------:|:-----------------:|
|    1     | 0 ≤ arg f < π/2   |
|    2     | π/2 ≤ arg f < π   |
|    3     | π ≤ arg f < 3π/2  |
|    4     | 3π/2 ≤ arg f < 2π |
"""
@inline function quadrant(val)
    # This function correponds to `vinq.m`
    rv, iv = reim(val)
    if (rv > 0) & (iv >= 0)
        return 1
    elseif (rv <= 0) & (iv > 0)
        return 2
    elseif (rv < 0) & (iv <= 0)
        return 3
    else
        # (rv >= 0) & (iv < 0)
        return 4
    end
end

"""
    assignquadrants!(quadrants, nodes, f, multithreading=false)

Evaluate function `f` for [`quadrant`](@ref) at `nodes` and fill `quadrants` in-place.

Each element of `quadrants` corresponds to the `index` of `IndexablePoint2D` in `nodes`.
"""
function assignquadrants!(quadrants, nodes, f::ScaledFunction{F,G}, multithreading=false) where {F,G}
    if multithreading
        @threads for p in nodes
            quadrants[getindex(p)] = quadrant(f(p))
        end
    else
        for p in nodes
            quadrants[getindex(p)] = quadrant(f(p))
        end
    end
    return nothing
end

"""
    candidateedges!(E, tess, quadrants)

Fill in `Vector` of candidate edges `E` that contain a phase change of 2 quadrants.

Any root or pole is located at the point where the regions described by four different
quadrants meet. Since any triangulation of the four nodes located in the four different
quadrants requires at least one edge of ``|ΔQ| = 2``, then all such edges are potentially in
the vicinity of a root or pole.

`E` is not sorted.
"""
function candidateedges!(E, tess, quadrants)
    # 10%+ better performance than even the v0.4.1 `delaunayedges()`
    # based on: https://github.com/JuliaGeometry/VoronoiDelaunay.jl/issues/47
    @inbounds for ix in 2:tess._last_trig_index
        tr = tess._trigs[ix]
        isexternal(tr) && continue

        # precalculate
        atr, btr, ctr = geta(tr), getb(tr), getc(tr)
        atri, btri, ctri = getindex(atr), getindex(btr), getindex(ctr)

        ix_na = tr._neighbour_a
        if ix_na > ix || isexternal(tess._trigs[ix_na])
            ΔQ = mod(quadrants[btri] - quadrants[ctri], 4)  # phase difference
            if ΔQ == 2
                edge = DelaunayEdge(btr, ctr)
                push!(E, edge)
            end
        end

        ix_nb = tr._neighbour_b
        if ix_nb > ix || isexternal(tess._trigs[ix_nb])
            ΔQ = mod(quadrants[atri] - quadrants[ctri], 4)
            if ΔQ == 2
                edge = DelaunayEdge(atr, ctr)
                push!(E, edge)
            end
        end

        ix_nc = tr._neighbour_c
        if ix_nc > ix || isexternal(tess._trigs[ix_nc])
            ΔQ = mod(quadrants[atri] - quadrants[btri], 4)
            if ΔQ == 2
                edge = DelaunayEdge(atr, btr)
                push!(E, edge)
            end
        end
    end

    return nothing
end

"""
    candidateedges!(E, phasediffs, tess, quadrants)

If `phasediffs` is a `Vector`, then the phase difference across each edge is `push!`ed
into `phasediffs`. This is useful for plotting.

Both `E` and `phasediffs` are updated in place.
"""
function candidateedges!(E, phasediffs, tess, quadrants)
    for edge in delaunayedges(tess)
        nodea, nodeb = geta(edge), getb(edge)
        idxa, idxb = getindex(nodea), getindex(nodeb)

        # NOTE: To match Matlab, force `idxa` < `idxb`
        # (order doesn't matter for `ΔQ == 2`, which is the only case we care about)
        # if idxa > idxb
        #     idxa, idxb = idxb, idxa
        # end

        ΔQ = mod(quadrants[idxa] - quadrants[idxb], 4)  # phase difference
        if ΔQ == 2
            push!(E, edge)
        end

        if phasediffs isa Vector
            push!(phasediffs, ΔQ)
        end
    end
end

"""
    zone(triangle, edge_idxs)

Return zone `1` or `2` for `DelaunayTriangle` `triangle`.

Zone `1` triangles have more than one node in `edge_idxs`, whereas zone `2` triangles have
only a single node.
"""
function zone(triangle, edge_idxs)
    na = geta(triangle)
    nb = getb(triangle)
    nc = getc(triangle)

    nai = getindex(na)
    nbi = getindex(nb)
    nci = getindex(nc)

    zone2 = false
    for idx in edge_idxs
        # with Julia 1.4.2: | is faster than || here
        if (nai == idx) | (nbi == idx) | (nci == idx)
            if !zone2
                # we need to keep searching b/c it might be zone 1
                zone2 = true
            else
                # zone1 = true
                return 1
            end
        end
    end

    return zone2 ? 2 : 0
end

"""
    uniqueindices!(idxs, edges)

Compute unique indices `idxs` in-place of all nodes in `edges`.
"""
function uniqueindices!(idxs, edges)
    empty!(idxs)
    for i in eachindex(edges)
        @inbounds ei = edges[i]
        push!(idxs, getindex(geta(ei)))
        push!(idxs, getindex(getb(ei)))
    end
    sort!(idxs)  # calling `sort!` first makes `unique!` more efficient
    unique!(idxs)
    return nothing
end

"""
    zone1newnodes!(newnodes, triangles, g2f, tolerance)

Add nodes (points) to `newnodes` in-place if they are in zone 1, i.e. triangles that had more than
one node.

`tolerance` is the minimum edge length an edge must have to go in `newnodes`.
"""
function zone1newnodes!(newnodes, triangles, g2f, tolerance)
    triangle1 = triangles[1]
    n1a = geta(triangle1)
    n1b = getb(triangle1)
    push!(newnodes, (n1a+n1b)/2)

    @inbounds for ii = 1:length(triangles)-1
        triangle = triangles[ii]
        na = geta(triangle)
        nb = getb(triangle)
        nc = getc(triangle)

        addnewnode!(newnodes, nb, nc, g2f, tolerance)
        addnewnode!(newnodes, nc, na, g2f, tolerance)
        addnewnode!(newnodes, geta(triangles[ii+1]), getb(triangles[ii+1]), g2f, tolerance)
    end
    te = triangles[end]
    na = geta(te)
    nb = getb(te)
    nc = getc(te)
    addnewnode!(newnodes, nb, nc, g2f, tolerance)
    addnewnode!(newnodes, nc, na, g2f, tolerance)

    # Remove the first of `newnodes` if the edge is too short
    distance(g2f(n1a), g2f(n1b)) < tolerance && popfirst!(newnodes)
    return nothing
end

@inline function addnewnode!(newnodes, node1, node2, g2f, tolerance)
    if distance(g2f(node1), g2f(node2)) > tolerance
        avgnode = (node1+node2)/2
        for p in newnodes
            distance(p, avgnode) < 2*eps() && return nothing
        end
        push!(newnodes, avgnode)  # only executed if we haven't already returned
    end
    return nothing
end

"""
    zone2newnodes!(newnodes, triangle, skinnytriangle)

Add node to `newnodes` for zone 2 ("skinny") triangles.

`skinnytriangle` is the maximum allowed ratio of the longest to shortest side length.
"""
@inline function zone2newnode!(newnodes, triangle, skinnytriangle)
    na = geta(triangle)
    nb = getb(triangle)
    nc = getc(triangle)

    # For skinny triangle check, `geom2fcn` not needed because units cancel out
    l1 = distance(na, nb)
    l2 = distance(nb, nc)
    l3 = distance(nc, na)
    if max(l1,l2,l3)/min(l1,l2,l3) > skinnytriangle
        avgnode = (na+nb+nc)/3
        push!(newnodes, avgnode)
    end
    return nothing
end

"""
    splittriangles!(zone1triangles, newnodes, tess, edge_idxs, params)

Add zone 2 triangles to `newnodes` and then update `Vector` `zone1triangles`, which require
special handling.
"""
function splittriangles!(zone1triangles, newnodes, tess, edge_idxs, params)
    empty!(zone1triangles)
    for triangle in tess
        z = zone(triangle, edge_idxs)

        if z == 1
            push!(zone1triangles, triangle)
        elseif z == 2
            zone2newnode!(newnodes, triangle, params.skinnytriangle)
        end
    end
    return nothing
end

"""
    findnextnode(prevnode, refnode, nodes, g2f)

Find the index of the next node in `nodes` as part of the candidate region boundary process.
The next one (after the reference) is picked from the fixed set of nodes.
"""
@inline function findnextnode(prevnode, refnode, nodes, g2f)
    P = g2f(prevnode)
    S = g2f(refnode)

    minphi = 2π + 1  # max diff of angles is 2π, so this is guaranteed larger
    minphi_idx = firstindex(nodes)

    for i in eachindex(nodes)
        N = g2f(nodes[i])

        SP = P - S
        SN = N - S

        phi = mod2pi(angle(SP) - angle(SN))

        if phi < minphi
            minphi = phi
            minphi_idx = i
        end
    end

    return minphi_idx
end

"""
    contouredges(tess, edges)

Find contour edges from all candidate `edges`.
"""
function contouredges(tess, edges)
    C = Vector{DelaunayEdge{IndexablePoint2D}}()
    sizehint!(C, length(edges))

    # Edges of triangles that contain at least 1 of `edges`
    for triangle in tess
        pa, pb, pc = geta(triangle), getb(triangle), getc(triangle)
        pai, pbi, pci = getindex(pa), getindex(pb), getindex(pc)

        for edge in edges
            eai, ebi = getindex(geta(edge)), getindex(getb(edge))

            if ((eai == pai) && (ebi == pbi)) | ((eai == pbi) && (ebi == pai)) |
                ((eai == pbi) && (ebi == pci)) | ((eai == pci) && (ebi == pbi)) |
                ((eai == pci) && (ebi == pai)) | ((eai == pai) && (ebi == pci))
                push!(C, DelaunayEdge(pa,pb), DelaunayEdge(pb,pc), DelaunayEdge(pc,pa))
                break  # only count each triangle once
            end
        end
    end

    # Remove duplicate edges
    sameunique!(C)

    return C
end

"""
    evaluateregions!(C, g2f)
"""
function evaluateregions!(C, g2f)
    # NOTE: The nodes of each region are in reverse order compared to Matlab with respect
    # to their quadrants

    # Initialize
    numregions = 1
    regions = [[geta(C[1])]]
    refnode = getb(C[1])
    popfirst!(C)

    nextedgeidxs = similar(Array{Int}, axes(C))
    empty!(nextedgeidxs)
    while length(C) > 0

        # This loop is equivalent to `findall(e->geta(e)==refnode, C)`
        # but avoids closure Core.Box issue
        for i in eachindex(C)
            if geta(C[i]) == refnode
                push!(nextedgeidxs, i)
            end
        end

        if !isempty(nextedgeidxs)
            if length(nextedgeidxs) == 1
                nextedgeidx = only(nextedgeidxs)
            else
                prevnode = regions[numregions][end]
                tempnodes = getb.(C[nextedgeidxs])
                idx = findnextnode(prevnode, refnode, tempnodes, g2f)
                nextedgeidx = nextedgeidxs[idx]
            end

            nextedge = C[nextedgeidx]
            push!(regions[numregions], geta(nextedge))
            refnode = getb(nextedge)
            deleteat!(C, nextedgeidx)
        else # isempty
            push!(regions[numregions], refnode)
            # New region
            numregions += 1
            push!(regions, [geta(C[1])])
            refnode = getb(C[1])
            popfirst!(C)
        end

        # reset `nextedgeidxs`
        empty!(nextedgeidxs)
    end

    push!(regions[numregions], refnode)

    return regions
end

"""
    rootsandpoles(regions, quadrants, g2f)

Identify roots and poles of function based on `regions` (usually a `Vector{Vector{IndexablePoint2D}}`)
and `quadrants`.
"""
function rootsandpoles(regions, quadrants, g2f::Geometry2Function{T}) where T
    complexT = complex(T)
    zroots = Vector{complexT}()
    zpoles = Vector{complexT}()
    for r in regions
        quadrantsequence = [quadrants[getindex(node)] for node in r]

        # Sign flip because `r` are in opposite order of Matlab?
        dquadrantsequence = -diff(quadrantsequence)
        for dq in dquadrantsequence
            if dq == 3
                dq = -1
            elseif dq == -3
                dq = 1
            elseif abs(dq) == 2
                # ``|ΔQ| = 2`` is ambiguous; cannot tell whether phase increases or
                # decreases by two quadrants
                dq = 0
            end
        end
        q = sum(dquadrantsequence)/4
        z = sum(g2f.(r))/length(r)

        if q > 0
            push!(zroots, convert(complexT, z))  # convert in case T isn't Float64
        elseif q < 0
            push!(zpoles, convert(complexT, z))
        end
    end

    return zroots, zpoles
end

"""
    tesselate!(tess, newnodes, f::ScaledFunction, params, pd=nothing)

Label quadrants, identify candidate edges, and iteratively split triangles, returning
the tuple `(tess, E, quadrants)`.
"""
function tesselate!(tess, newnodes, f::ScaledFunction, params, pd=nothing)
    # Initialize
    numnodes = tess._total_points_added
    @assert numnodes == 0

    g2f = f.g2f

    E = Vector{DelaunayEdge{IndexablePoint2D}}()
    phasediffs = Vector{Int}()
    quadrants = Vector{Int}()
    edge_idxs = Vector{Int}()
    zone1triangles = Vector{DelaunayTriangle{IndexablePoint2D}}()

    iteration = 0
    while (iteration < params.maxiterations) && (numnodes < params.maxnodes)
        iteration += 1

        # Determine which quadrant function value belongs at each node
        numnewnodes = length(newnodes)
        append!(quadrants, Vector{Int}(undef, numnewnodes))
        assignquadrants!(quadrants, newnodes, f, params.multithreading)

        # Add new nodes to `tess`
        push!(tess, newnodes)
        numnodes += numnewnodes

        # Determine candidate edges that may be near a root or pole
        empty!(E)  # start with a blank E
        if pd isa PlotData
            empty!(phasediffs)
            candidateedges!(E, phasediffs, tess, quadrants)
        else
            candidateedges!(E, tess, quadrants)
        end
        isempty(E) && return tess, E, quadrants, phasediffs  # no roots or poles found

        # Select candidate edges that are longer than the chosen tolerance
        selectE = filter(e -> longedge(e, params.tolerance, g2f), E)
        isempty(selectE) && return tess, E, quadrants, phasediffs

        # return if maximum edge length has reached tolerance
        maxElength = maximum(distance(g2f(e)) for e in selectE)
        maxElength < params.tolerance && return tess, E, quadrants, phasediffs

        # Get unique indices of nodes in `edges`
        uniqueindices!(edge_idxs, selectE)

        # Refine (split) triangles
        empty!(newnodes)
        splittriangles!(zone1triangles, newnodes, tess, edge_idxs, params)

        # Add new nodes in zone 1
        zone1newnodes!(newnodes, zone1triangles, g2f, params.tolerance)

        # Have to assign indexes to new nodes (which are all currently -1)
        setindex!.(newnodes, (1:length(newnodes)).+numnodes)
    end

    (iteration >= params.maxiterations) && @warn "params.maxiterations reached"
    (numnodes >= params.maxnodes) && @warn "params.maxnodes reached"

    return tess, E, quadrants, phasediffs
end

function setup(origcoords)
    rmin, rmax = minimum(real, origcoords), maximum(real, origcoords)
    imin, imax = minimum(imag, origcoords), maximum(imag, origcoords)

    # Scaling parameters
    g2f = Geometry2Function(rmin, rmax, imin, imax)
    scaledorigcoords = fcn2geom.(origcoords, rmin, rmax, imin, imax)

    @assert minimum(real, scaledorigcoords) >= min_coord &&
        minimum(imag, scaledorigcoords) >= min_coord &&
        maximum(real, scaledorigcoords) <= max_coord &&
        maximum(imag, scaledorigcoords) <= max_coord "Scaled coordinates out of bounds"

    newnodes = [IndexablePoint2D(real(coord), imag(coord), idx) for (idx, coord) in
                enumerate(scaledorigcoords)]

    return newnodes, g2f
end

"""
    grpf(fcn, origcoords, params=GRPFParams())

Return a vector `roots` and a vector `poles` of a single (complex) argument function
`fcn`.

Searches within a domain specified by the vector of complex `origcoords`.

# Examples
```jldoctest
julia> simplefcn(z) = (z - 1)*(z - im)^2*(z + 1)^3/(z + im)

julia> xb, xe = -2, 2

julia> yb, ye = -2, 2

julia> r = 0.1

julia> origcoords = rectangulardomain(complex(xb, yb), complex(xe, ye), r)

julia> roots, poles = grpf(simplefcn, origcoords);

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
function grpf(fcn, origcoords, params=GRPFParams())
    newnodes, g2f = setup(origcoords)
    f = ScaledFunction(fcn, g2f)

    tess = DelaunayTessellation2D{IndexablePoint2D}(params.tess_sizehint)

    tess, E, quadrants, _ = tesselate!(tess, newnodes, f, params)

    complexT = complex(eltype(g2f))
    isempty(E) && return Vector{complexT}(), Vector{complexT}()

    C = contouredges(tess, E)
    regions = evaluateregions!(C, g2f)
    zroots, zpoles = rootsandpoles(regions, quadrants, g2f)

    return zroots, zpoles
end

"""
    grpf(fcn, origcoords, ::PlotData, params=GRPFParams())

Variant of `grpf` that returns `quadrants`, `phasediffs`, the `VoronoiDelauany`
tesselation `tess`, and the `Geometry2Function` struct `g2f` for converting from the
`VoronoiDelaunay` to function space, in addition to `zroots` and `zpoles`.

These additional outputs are primarily for plotting or diagnostics.

# Examples
```
julia> roots, poles, quadrants, phasediffs, tess, g2f = grpf(simplefcn, origcoords, PlotData());
```
"""
function grpf(fcn, origcoords, ::PlotData, params=GRPFParams())
    newnodes, g2f = setup(origcoords)
    f = ScaledFunction(fcn, g2f)

    tess = DelaunayTessellation2D{IndexablePoint2D}(params.tess_sizehint)

    tess, E, quadrants, phasediffs = tesselate!(tess, newnodes, f, params, PlotData())

    complexT = complex(eltype(g2f))
    isempty(E) && return (Vector{complexT}(), Vector{complexT}(), quadrants, phasediffs,
                          tess, g2f)

    C = contouredges(tess, E)
    regions = evaluateregions!(C, g2f)
    zroots, zpoles = rootsandpoles(regions, quadrants, g2f)

    return zroots, zpoles, quadrants, phasediffs, tess, g2f
end

end # module
