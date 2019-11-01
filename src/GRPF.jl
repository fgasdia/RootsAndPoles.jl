__precompile__()

"""
# GRPF: Global complex Roots and Poles Finding algorithm

A Julia implementation of the GRPF (https://github.com/PioKow/GRPF) by Piotr
Kowalczyk. Matlab files ('.m') listed in _see also_ sections of function doc strings refer
to this git repo.
"""
module GRPF

using LinearAlgebra
using StaticArrays
import GeometricalPredicates: intriangle
using VoronoiDelaunay

include("VoronoiDelaunayExtensions.jl")
include("GeneralFunctions.jl")

export rectangulardomain, diskdomain, grpf, grpf!, PhaseDiffs

const MAXITERATIONS = 100
const MAXNODES = 500000
const SKINNYTRIANGLE = 3

struct PhaseDiffs end

"""
    quadrant(val)

Convert complex function value `val` to quadrant number.

See also: `vinq.m`
"""
@inline function quadrant(val::Complex)::Int8
    rv, iv = reim(val)
    if (rv > 0) & (iv >= 0)
        return 1
    elseif (rv <= 0) & (iv > 0)
        return 2
    elseif (rv < 0) & (iv <= 0)
        return 3
    elseif (rv >= 0) & (iv < 0)
        return 4
    else
        error("Function value $val cannot be assigned to quadrant.")
    end
end

"""
    assignquadrants!(quadrants, nodes, fcn)

Evaluate `fcn` for [`quadrant`](@ref) at `nodes` and fill `quadrants`.

`quadrants` is a Vector{} where each index corresponds to `node` index.
"""
@inline function assignquadrants!(quadrants::Vector{Int8},
                          nodes::Vector{IndexablePoint2D}, fcn::Function)
    @inbounds for ii in eachindex(nodes)
        val = fcn(nodes[ii])  # TODO: `val` is of Any type
        quadrants[getindex(nodes[ii])] = quadrant(complex(val))
    end
    nothing
end

"""
    candidateedges(tess, quadrants)

Return candidate edges `𝓔` that contain a phase change of 2 quadrants.

Any root or pole is located at the point where the regions described by four different
quadrants meet. Since any triangulation of the four nodes located in the four different
quadrants requires at least one edge of ``|ΔQ| = 2``, then all such edges are potentially
in the vicinity of a root or pole.

Notes:
 - Order of `𝓔` is not guaranteed.
 - Count of phasediffs `ΔQ` of value 1 and 3 can differ from Matlab in normal operation,
 because it depends on "direction" of edge.
"""
function candidateedges(
    tess::DelaunayTessellation2D{IndexablePoint2D},
    quadrants::Vector{Int8}
    )

    𝓔 = Vector{DelaunayEdge{IndexablePoint2D}}()

    @inbounds for edge in delaunayedges_fast(tess)
        e = edge::DelaunayEdge{IndexablePoint2D}
        nodea, nodeb = geta(e), getb(e)
        idxa, idxb = getindex(nodea), getindex(nodeb)

        # NOTE: To match Matlab, force `idxa` < `idxb`
        # (order doesn't matter for `ΔQ == 2`, which is the only case we care about)
        # if idxa > idxb
        #     idxa, idxb = idxb, idxa
        # end

        ΔQ = mod(quadrants[idxa] - quadrants[idxb], 4)  # phase difference
        if ΔQ == 2
            push!(𝓔, e)
        end
    end
    return 𝓔
end

function candidateedges(
    tess::DelaunayTessellation2D{IndexablePoint2D},
    quadrants::Vector{Int8},
    ::PhaseDiffs
    )

    𝓔 = Vector{DelaunayEdge{IndexablePoint2D}}()
    phasediffs = Vector{Int8}()

    @inbounds for edge in delaunayedges_fast(tess)
        e = edge::DelaunayEdge{IndexablePoint2D}
        nodea, nodeb = geta(e), getb(e)
        idxa, idxb = getindex(nodea), getindex(nodeb)

        # NOTE: To match Matlab, force `idxa` < `idxb`
        # (order doesn't matter for `ΔQ == 2`, which is the only case we care about)
        if idxa > idxb
            idxa, idxb = idxb, idxa
        end

        ΔQ = mod(quadrants[idxa] - quadrants[idxb], 4)  # phase difference
        if ΔQ == 2
            push!(𝓔, e)
        end

        push!(phasediffs, ΔQ)
    end
    return 𝓔, phasediffs
end

"""
    counttriangleswithnodes(tess, edges)

Count how many times each triangle contains a node in `edges`.
"""
function counttriangleswithnodes(
    tess::DelaunayTessellation2D{IndexablePoint2D},
    edges::Vector{DelaunayEdge{IndexablePoint2D}}
    )

    # Nodes of select edges
    edgenodes = Vector{IndexablePoint2D}()
    uniquenodes!(edgenodes, edges)

    trianglecounts = zeros(Int, count(.!isexternal.(tess._trigs)))
    triidx = 0
    @inbounds for triangle in tess
        triidx += 1
        # `triangle` is in general not equal to `tess._trigs[triidx]`
        ea = geta(triangle)
        eb = getb(triangle)
        ec = getc(triangle)
        @inbounds for nodeidx in eachindex(edgenodes)
            if (ea == edgenodes[nodeidx]) || (eb == edgenodes[nodeidx]) || (ec == edgenodes[nodeidx])
                trianglecounts[triidx] += 1
            end
        end
    end
    return trianglecounts
end

function uniquenodes!(
    edgenodes::Vector{IndexablePoint2D},
    edges::Vector{DelaunayEdge{IndexablePoint2D}}
    )

    @inbounds for ii in eachindex(edges)
        nodea = geta(edges[ii])
        nodeb = getb(edges[ii])

        nodea in edgenodes || push!(edgenodes, nodea)
        nodeb in edgenodes || push!(edgenodes, nodeb)
    end
    nothing
end

"""
    zone1newnodes!(newnodes, triangles, geom2fcn, tolerance)

Add nodes to `newnodes` in zone 1, i.e. triangles that had more than one node.
"""
function zone1newnodes!(
    newnodes::Vector{IndexablePoint2D},
    triangles::Vector{DelaunayTriangle{IndexablePoint2D}},
    geom2fcn::Function,
    tolerance
    )

    triangle1 = triangles[1]
    n1a = geta(triangle1)
    n1b = getb(triangle1)
    push!(newnodes, (n1a+n1b)/2)

    @inbounds for ii = 1:length(triangles)-1
        na = geta(triangles[ii])
        nb = getb(triangles[ii])
        nc = getc(triangles[ii])

        addnewnode!(newnodes, nb, nc, geom2fcn, tolerance)
        addnewnode!(newnodes, nc, na, geom2fcn, tolerance)
        addnewnode!(newnodes, geta(triangles[ii+1]), getb(triangles[ii+1]), geom2fcn, tolerance)
    end
    na = geta(triangles[end])
    nb = getb(triangles[end])
    nc = getc(triangles[end])
    addnewnode!(newnodes, nb, nc, geom2fcn, tolerance)
    addnewnode!(newnodes, nc, na, geom2fcn, tolerance)

    # Remove the first of `newnodes` if the edge is too short
    distance(geom2fcn(n1a), geom2fcn(n1b)) < tolerance && popfirst!(newnodes)
    nothing
end

@inline function addnewnode!(
    newnodes::Vector{IndexablePoint2D},
    node1::IndexablePoint2D,
    node2::IndexablePoint2D,
    geom2fcn::Function,
    tolerance
    )

    if distance(geom2fcn(node1), geom2fcn(node2)) > tolerance
        avgnode = (node1+node2)/2
        @inbounds for ii in eachindex(newnodes)
            distance(newnodes[ii], avgnode) < 2*eps() && return nothing
        end
        push!(newnodes, avgnode)  # only executed if we haven't already returned
    end
    nothing
end

"""
    zone2newnodes!(newnodes, triangles)

Add nodes to `newnodes` in zone 2 (skinny triangles).
"""
@inline function zone2newnodes!(
    newnodes::Vector{IndexablePoint2D},
    triangles::Vector{DelaunayTriangle{IndexablePoint2D}}
    )

    @inbounds for triangle in triangles
        na = geta(triangle)
        nb = getb(triangle)
        nc = getc(triangle)

        # For skinny triangle check, `geom2fcn` not needed because units cancel out
        l1 = distance(na, nb)
        l2 = distance(nb, nc)
        l3 = distance(nc, na)
        if max(l1,l2,l3)/min(l1,l2,l3) > SKINNYTRIANGLE
            avgnode = (na+nb+nc)/3
            push!(newnodes, avgnode)
        end
    end
    nothing
end

"""
    contouredges(tess, edges)

Find contour edges from all candidate edges.
"""
function contouredges(
    tess::DelaunayTessellation2D{IndexablePoint2D},
    edges::Vector{DelaunayEdge{IndexablePoint2D}}
    )

    # Edges of triangles that contain at least 1 of `edges`
    tmpedges = Vector{DelaunayEdge{IndexablePoint2D}}()
    @inbounds for triangle in tess
        # We don't know which "direction" the edges are defined in the triangle,
        # so we need to test both
        pa, pb, pc = geta(triangle), getb(triangle), getc(triangle)

        edgea = DelaunayEdge(pa, pb)
        edgearev = DelaunayEdge(pb, pa)
        edgeb = DelaunayEdge(pb, pc)
        edgebrev = DelaunayEdge(pc, pb)
        edgec = DelaunayEdge(pc, pa)
        edgecrev = DelaunayEdge(pa, pc)

        # Does triangle contain edge?
        @inbounds for edge in edges
            if edgea == edge || edgeb == edge || edgec == edge ||
                edgearev == edge || edgebrev == edge || edgecrev == edge
                push!(tmpedges, edgea, edgeb, edgec)
                break  # only count each triangle once
            end
        end
    end

    # Remove duplicate (reverse) edges from `tmpedges` and otherwise append to `𝐶`
    𝐶 = Vector{DelaunayEdge{IndexablePoint2D}}()
    duplicateedges = zeros(Int8, length(tmpedges))
    @inbounds for (idxa, edgea) in enumerate(tmpedges)
        if duplicateedges[idxa] == 0
            @inbounds for (idxb, edgeb) in enumerate(tmpedges)
                # Check if Edge(a,b) == Edge(b, a), i.e. there are duplicate edges
                if edgea == DelaunayEdge(getb(edgeb), geta(edgeb))
                    duplicateedges[idxa] = 2
                    duplicateedges[idxb] = 2
                    break
                end
            end
            if duplicateedges[idxa] != 2
                duplicateedges[idxa] = 1
                push!(𝐶, edgea)
            end
        end
    end

    return 𝐶
end

"""
    splittriangles(tess, trianglecounts)

Separate triangles in `tess` by zones.
"""
function splittriangles(
    tess::DelaunayTessellation2D{IndexablePoint2D},
    trianglecounts::Vector{Int}
    )

    zone1triangles = Vector{DelaunayTriangle{IndexablePoint2D}}()
    zone2triangles = Vector{DelaunayTriangle{IndexablePoint2D}}()
    ii = 0
    @inbounds for triangle in tess
        ii += 1
        if trianglecounts[ii] > 1
            push!(zone1triangles, triangle)
        elseif trianglecounts[ii] == 1
            push!(zone2triangles, triangle)
        end
    end
    return zone1triangles, zone2triangles
end

"""
    findnextnode(prevnode, refnode, tempnodes, geom2fcn)

Find the index of the next node in the candidate region boundary process. The next one (after
the reference) is picked from the fixed set of nodes.

See also: `FindNextNode.m`
"""
function findnextnode(
    prevnode::IndexablePoint2D,
    refnode::IndexablePoint2D,
    tempnodes::Vector{IndexablePoint2D},
    geom2fcn::Function
    )

    P = geom2fcn(prevnode)
    S = geom2fcn(refnode)

    ϕs = Vector{Float64}(undef, length(tempnodes))
    for i in eachindex(tempnodes)
        N = geom2fcn(tempnodes[i])

        SP = P - S
        SN = N - S

        SPlength = norm(SP)
        SNlength = norm(SN)

        dotprod = real(SP)*real(SN) + imag(SP)*imag(SN)
        ϕs[i] = acos(dotprod/(SPlength*SNlength))
        if real(SP)*imag(SN) - imag(SP)*real(SN) < 0
            ϕs[i] = 2π - ϕs[i]
        end
    end
    return findmin(ϕs)[2]  # return index of minimum `ϕ`
end

"""
    evaluateregions!(𝐶, geom2fcn)

TODO: Go in reverse (from `end` rather than `1` so we don't popfirst!) or use careful indexing
and don't pop at all

# Note: The nodes of each region are in reverse order compared to Matlab wrt their quadrants?
"""
function evaluateregions!(
    𝐶::Vector{DelaunayEdge{IndexablePoint2D}},
    geom2fcn::Function
    )

    # Initialize
    numregions = 1

    # TODO: Type instability of `refnode` > Core.Box
    # see https://github.com/JuliaLang/julia/issues/15276#issuecomment-297596373

    regions = [[geta(𝐶[1])]]
    refnode = getb(𝐶[1])
    popfirst!(𝐶)
    while length(𝐶) > 0
        nextedgeidxs = findall(e->geta(e)==refnode, 𝐶)
        if isempty(nextedgeidxs)
            push!(regions[numregions], refnode)
            # New region
            numregions += 1
            push!(regions, [geta(𝐶[1])])
            refnode = getb(𝐶[1])
            popfirst!(𝐶)
        else
            if length(nextedgeidxs) > 1
                prevnode = regions[numregions][end]
                tempnodes = getb.(𝐶[nextedgeidxs])
                idx = findnextnode(prevnode, refnode, tempnodes, geom2fcn)
                nextedgeidx = nextedgeidxs[idx]
            else
                nextedgeidx = nextedgeidxs[1]
            end

            nextedge = 𝐶[nextedgeidx]
            push!(regions[numregions], geta(nextedge))
            refnode = getb(nextedge)
            deleteat!(𝐶, nextedgeidx)
        end
    end
    push!(regions[numregions], refnode)

    return regions
end

"""
    rootsandpoles(regions, quadrants, geom2fcn)

Identify roots and poles of function based on regions and quadrants.
"""
function rootsandpoles(
    regions::Vector{Vector{IndexablePoint2D}},
    quadrants::Vector{Int8},
    geom2fcn::Function
    )

    numregions = length(regions)

    zroots = Vector{ComplexF64}()
    zpoles = Vector{ComplexF64}()
    @inbounds for ii in eachindex(regions)
        # XXX: Order of regions?
        quadrantsequence = [quadrants[getindex(node)] for node in regions[ii]]

        # Sign flip because `regions[ii]` are in opposite order of Matlab?
        dQ = -diff(quadrantsequence)
        @inbounds for jj in eachindex(dQ)
            if dQ[jj] == 3
                dQ[jj] = -1
            elseif dQ[jj] == -3
                dQ[jj] = 1
            elseif abs(dQ[jj]) == 2
                # ``|ΔQ| = 2`` is ambiguous; cannot tell whether phase increases or decreases by two quadrants
                dQ[jj] = 0
            end
        end
        q = sum(dQ)/4
        z = sum(geom2fcn.(regions[ii]))/length(regions[ii])

        if q > 0
            push!(zroots, z)
        elseif q < 0
            push!(zpoles, z)
        end
    end

    return zroots, zpoles
end

"""
    tesselate!(tess, newnodes, fcn, geom2fcn, tolerance)

Label quadrants, identify candidate edges, and iteratively split triangles.
"""
function tesselate!(
    tess::DelaunayTessellation2D{IndexablePoint2D},
    newnodes::Vector{IndexablePoint2D},
    fcn::Function,
    geom2fcn::Function,
    tolerance
    )

    # Initialize
    numnodes = tess._total_points_added
    @assert numnodes == 0

    𝓔 = Vector{DelaunayEdge{IndexablePoint2D}}()
    quadrants = Vector{Int8}()

    iteration = 0
    while (iteration < MAXITERATIONS) && (numnodes < MAXNODES)
        iteration += 1

        # Determine which quadrant function value belongs at each node
        numnewnodes = length(newnodes)
        append!(quadrants, Vector{Int8}(undef, numnewnodes))
        assignquadrants!(quadrants, newnodes, fcn)

        # Add new nodes to `tess`
        push!(tess, newnodes)
        numnodes += numnewnodes

        # Determine candidate edges that may be near a root or pole
        𝓔 = candidateedges(tess, quadrants)
        isempty(𝓔) && error("No roots in the domain")

        # Select candidate edges that are longer than the chosen tolerance
        select𝓔 = filter(e -> longedge(e, tolerance, geom2fcn), 𝓔)
        isempty(select𝓔) && return tess, 𝓔, quadrants

        max𝓔length = maximum(distance(geom2fcn(e)) for e in select𝓔)
        max𝓔length < tolerance && return tess, 𝓔, quadrants, phasediffs

        # How many times does each triangle contain a `select𝓔` node?
        trianglecounts = counttriangleswithnodes(tess, select𝓔)
        zone1triangles, zone2triangles = splittriangles(tess, trianglecounts)

        # Add new nodes in zone 1
        newnodes = Vector{IndexablePoint2D}()
        zone1newnodes!(newnodes, zone1triangles, geom2fcn, tolerance)

        # Add new nodes in zone 2
        zone2newnodes!(newnodes, zone2triangles)

        # Have to assign indexes to new nodes (which are all currently -1)
        setindex!.(newnodes, (1:length(newnodes)).+numnodes)
    end

    return tess, 𝓔, quadrants
end

function tesselate!(
    tess::DelaunayTessellation2D{IndexablePoint2D},
    newnodes::Vector{IndexablePoint2D},
    fcn::Function,
    geom2fcn::Function,
    tolerance,
    ::PhaseDiffs
    )

    # Initialize
    numnodes = tess._total_points_added
    @assert numnodes == 0

    𝓔 = Vector{DelaunayEdge{IndexablePoint2D}}()
    quadrants = Vector{Int8}()

    iteration = 0
    while (iteration < MAXITERATIONS) && (numnodes < MAXNODES)
        iteration += 1

        # Determine which quadrant function value belongs at each node
        numnewnodes = length(newnodes)
        append!(quadrants, Vector{Int8}(undef, numnewnodes))
        assignquadrants!(quadrants, newnodes, fcn)

        # Add new nodes to `tess`
        push!(tess, newnodes)
        numnodes += numnewnodes

        # Determine candidate edges that may be near a root or pole
        𝓔, phasediffs = candidateedges(tess, quadrants, PhaseDiffs())
        isempty(𝓔) && error("No roots in the domain")

        # Select candidate edges that are longer than the chosen tolerance
        select𝓔 = filter(e -> longedge(e, tolerance, geom2fcn), 𝓔)
        isempty(select𝓔) && return tess, 𝓔, quadrants, phasediffs

        max𝓔length = maximum(distance(geom2fcn(e)) for e in select𝓔)
        max𝓔length < tolerance && return tess, 𝓔, quadrants, phasediffs

        # How many times does each triangle contain a `select𝓔` node?
        trianglecounts = counttriangleswithnodes(tess, select𝓔)
        zone1triangles, zone2triangles = splittriangles(tess, trianglecounts)

        # Add new nodes in zone 1
        newnodes = Vector{IndexablePoint2D}()
        zone1newnodes!(newnodes, zone1triangles, geom2fcn, tolerance)

        # Add new nodes in zone 2
        zone2newnodes!(newnodes, zone2triangles)

        # Have to assign indexes to new nodes (which are all currently -1)
        setindex!.(newnodes, (1:length(newnodes)).+numnodes)
    end

    return tess, 𝓔, quadrants, phasediffs
end

"""
    grpf(fcn, origcoords, tolerance, tess_size_hint=5000)

Return roots and poles of a single (complex) argument function `fcn`.

Searches within a domain specified by the vector of `origcoords` with a final `tolerance`.
`tess_size_hint` is a sizehint for the DelaunayTessellation.

# Examples
```jldoctest
julia> simplefcn(z) = (z - 1)*(z - im)^2*(z + 1)^3/(z + im)

julia> xb, xe = -2, 2

julia> yb, ye = -2, 2

julia> r = 0.1

julia> tolerance = 1e-9

julia> origcoords = rectangulardomain(complex(xb, yb), complex(xe, ye), r)

julia> roots, poles = grpf(simplefcn, origcoords, tolerance);

julia> roots
3-element Array{Complex{Float64},1}:
    -0.9999999999512241 - 2.865605189037104e-11im
     0.9999999996829548 - 6.208811242913729e-11im
 1.9022756703179778e-10 + 1.0000000000372526im

 julia> poles
 1-element Array{Complex{Float64},1}:
 -3.8045513406359555e-10 - 1.0000000002235174im
```

See also: [`grpf`](@ref)
"""
function grpf(fcn::Function, origcoords::AbstractArray, tolerance, tess_size_hint=5000)
    # Need to map space domain for VoronoiDelaunay.jl
    rmin, rmax = minimum(real(origcoords)), maximum(real(origcoords))
    imin, imax = minimum(imag(origcoords)), maximum(imag(origcoords))

    # `max_coord` and `min_coord` are provided by `VoronoiDelaunay.jl`
    width = max_coord - min_coord
    ra = width/(rmax-rmin)
    rb = max_coord - ra*rmax

    ia = width/(imax-imin)
    ib = max_coord - ia*imax

    origcoords = mapfunctionval.(origcoords, ra, rb, ia, ib)
    newnodes = [IndexablePoint2D(real(coord), imag(coord), idx) for (idx, coord) in enumerate(origcoords)]
    tess = DelaunayTessellation2D{IndexablePoint2D}(tess_size_hint)

    f = pt -> fcn(geom2fcn(pt, ra, rb, ia, ib))
    g = z -> geom2fcn(z, ra, rb, ia, ib)

    tess, 𝓔, quadrants = tesselate!(tess, newnodes, f, g, tolerance)
    𝐶 = contouredges(tess, 𝓔)
    regions = evaluateregions!(𝐶, g)
    zroots, zpoles = rootsandpoles(regions, quadrants, g)

    return zroots, zpoles
end

"""
    grpf(fcn, origcoords, tolerance, ::PhaseDiffs, tess_size_hint=5000)

Variant of `grpf` that returns `quadrants` and `phasediffs` in addition to `zroots` and
`zpoles`, primarily for plotting or diagnostics.

# Examples
```jldoctest
julia> simplefcn(z) = (z - 1)*(z - im)^2*(z + 1)^3/(z + im)

julia> xb, xe = -2, 2

julia> yb, ye = -2, 2

julia> r = 0.1

julia> tolerance = 1e-9

julia> origcoords = rectangulardomain(complex(xb, yb), complex(xe, ye), r)

julia> roots, poles, quadrants, phasediffs = grpf(simplefcn, origcoords, tolerance, PhaseDiffs());

julia> roots
3-element Array{Complex{Float64},1}:
    -0.9999999999512241 - 2.865605189037104e-11im
     0.9999999996829548 - 6.208811242913729e-11im
 1.9022756703179778e-10 + 1.0000000000372526im

 julia> poles
 1-element Array{Complex{Float64},1}:
 -3.8045513406359555e-10 - 1.0000000002235174im
```

See also: [`grpf`](@ref)
"""
function grpf(fcn::Function, origcoords::AbstractArray, tolerance, ::PhaseDiffs, tess_size_hint=5000)
    # Need to map space domain for VoronoiDelaunay.jl
    rmin, rmax = minimum(real(origcoords)), maximum(real(origcoords))
    imin, imax = minimum(imag(origcoords)), maximum(imag(origcoords))

    # `max_coord` and `min_coord` are provided by `VoronoiDelaunay.jl`
    width = max_coord - min_coord
    ra = width/(rmax-rmin)
    rb = max_coord - ra*rmax

    ia = width/(imax-imin)
    ib = max_coord - ia*imax

    origcoords = mapfunctionval.(origcoords, ra, rb, ia, ib)
    newnodes = [IndexablePoint2D(real(coord), imag(coord), idx) for (idx, coord) in enumerate(origcoords)]
    tess = DelaunayTessellation2D{IndexablePoint2D}(tess_size_hint)

    f = pt -> fcn(geom2fcn(pt, ra, rb, ia, ib))
    g = z -> geom2fcn(z, ra, rb, ia, ib)

    tess, 𝓔, quadrants, phasediffs = tesselate!(tess, newnodes, f, g, tolerance, PhaseDiffs())
    𝐶 = contouredges(tess, 𝓔)
    regions = evaluateregions!(𝐶, g)
    zroots, zpoles = rootsandpoles(regions, quadrants, g)

    return zroots, zpoles, quadrants, phasediffs
end

end # module
