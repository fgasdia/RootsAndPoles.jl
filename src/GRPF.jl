__precompile__()

"""
# GRPF: Global complex Roots and Poles Finding algorithm

A Julia implementation of the GRPF (https://github.com/PioKow/GRPF) by Piotr
Kowalczyk. Matlab files ('.m') listed in _see also_ sections of function doc strings refer
to this git repo.
"""
module GRPF

using LinearAlgebra
using Statistics
using StaticArrays
import GeometricalPredicates: intriangle
using VoronoiDelaunay

include("VoronoiDelaunayExtensions.jl")
include("GeneralFunctions.jl")

export rectangulardomain, diskdomain, grpf

const maxiterations = 100
const maxnodes = 500000
const skinnytriangle = 3

"""
    quadrant(val)

Convert complex function value `val` to quadrant number.

See also: `vinq.m`
"""
function quadrant(val::Complex)::UInt8
    if (real(val) > 0) & (imag(val) >= 0)
        return 1
    elseif (real(val) <= 0) & (imag(val) > 0)
        return 2
    elseif (real(val) < 0) & (imag(val) <= 0)
        return 3
    elseif (real(val) >= 0) & (imag(val) < 0)
        return 4
    else
        # val = 0?
        error("Function value $val cannot be assigned to quadrant.")
    end
end

"""
    assignquadrants!(quadrants, nodes, fcn)

Evaluate `fcn` for [`quadrant`](@ref) at `nodes` and fill `quadrants`.

`quadrants` is a Vector{} where each index corresponds to `node` index.
"""
function assignquadrants!(quadrants::Vector{UInt8},
                          nodes::Vector{IndexablePoint2D}, fcn::Function)
    for ii in eachindex(nodes)
        val = complex(fcn(nodes[ii]))
        quadrants[getindex(nodes[ii])] = quadrant(val)
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
 - Count of `phasediffs` of value 1 and 3 can differ from Matlab in normal operation,
 because it depends on "direction" of edge.
 - `phasediffs` is only needed for diagnosis and plotting.
"""
function candidateedges(
    tess::DelaunayTessellation2D{IndexablePoint2D},
    quadrants::Vector{UInt8}
    )

    phasediffs = Vector{UInt8}()
    𝓔 = Vector{DelaunayEdge{IndexablePoint2D}}()

    edgeiter = delaunayedges(tess)
    for edge in edgeiter
        e::DelaunayEdge{IndexablePoint2D} = edge  # TODO: infer edge type
        nodea, nodeb = geta(e), getb(e)
        idxa, idxb = getindex(nodea), getindex(nodeb)

        # To match Matlab, force `idxa` < `idxb`
        # (order doesn't matter for `ΔQ == 2`, which is the only case we care about)
        if idxa > idxb
            idxa, idxb = idxb, idxa
        end

        ΔQ = mod(quadrants[idxa] - quadrants[idxb], 4)
        ΔQ == 2 && push!(𝓔, e)

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
    for triangle in tess
        triidx += 1
        # `triangle` is in general not equal to `tess._trigs[triidx]`
        ea = geta(triangle)
        eb = getb(triangle)
        ec = getc(triangle)
        for nodeidx in eachindex(edgenodes)
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

    for ii in eachindex(edges)
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
    for ii = 1:length(triangles)-1
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

function addnewnode!(
    newnodes::Vector{IndexablePoint2D},
    node1::IndexablePoint2D,
    node2::IndexablePoint2D,
    geom2fcn::Function,
    tolerance
    )

    if distance(geom2fcn(node1), geom2fcn(node2)) > tolerance
        avgnode = (node1+node2)/2
        for ii in eachindex(newnodes)
            distance(newnodes[ii], avgnode) < 2eps() && return nothing
        end
        push!(newnodes, avgnode)  # only executed if we haven't already returned
    end
    nothing
end

"""
Add nodes to `newnodes` in zone 2 (skinny triangles).
"""
function zone2newnodes!(
    newnodes::Vector{IndexablePoint2D},
    triangles::Vector{DelaunayTriangle{IndexablePoint2D}}
    )

    for ii in eachindex(triangles)
        triangle = triangles[ii]::DelaunayTriangle{IndexablePoint2D}
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
    for triangle in tess
        # We don't know which "direction" the edges are defined in the triangle,
        # so we need to test both
        edgea = DelaunayEdge(geta(triangle), getb(triangle))
        edgearev = DelaunayEdge(getb(triangle), geta(triangle))
        edgeb = DelaunayEdge(getb(triangle), getc(triangle))
        edgebrev = DelaunayEdge(getc(triangle), getb(triangle))
        edgec = DelaunayEdge(getc(triangle), geta(triangle))
        edgecrev = DelaunayEdge(geta(triangle), getc(triangle))

        # Does triangle contain edge?
        for edge in edges
            if edgea == edge || edgeb == edge || edgec == edge ||
                edgearev == edge || edgebrev == edge || edgecrev == edge
                push!(tmpedges, edgea, edgeb, edgec)
                break  # only count each triangle once
            end
        end
    end

    # Remove duplicate (reverse) edges from `tmpedges` and otherwise append to `𝐶`
    𝐶 = Vector{DelaunayEdge{IndexablePoint2D}}()
    duplicateedges = zeros(Int, length(tmpedges))
    for (idxa, edgea) in enumerate(tmpedges)
        if duplicateedges[idxa] == 0
            for (idxb, edgeb) in enumerate(tmpedges)
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

Separate triangles by zones.
"""
function splittriangles(
    tess::DelaunayTessellation2D{IndexablePoint2D},
    trianglecounts::Vector{Int}
    )

    zone1triangles = Vector{DelaunayTriangle{IndexablePoint2D}}()
    zone2triangles = Vector{DelaunayTriangle{IndexablePoint2D}}()
    ii = 0
    for triangle in tess
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
    N = geom2fcn.(tempnodes)

    numtempnodes = length(N)

    SP = ones(numtempnodes)*(P-S)
    SN = N .- S

    SPlength = norm.(SP)
    SNlength = norm.(SN)

    dotprod = real.(SP).*real.(SN)+imag.(SP).*imag.(SN)
    ϕ = acos.(dotprod./(SPlength.*SNlength))
    tmp = findall(real.(SP).*imag.(SN)-imag.(SP).*real.(SN) .< 0)
    ϕ[tmp] .= 2π .- ϕ[tmp]
    findmin(ϕ)[2]  # return index of minimum `ϕ`
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

    regions = [[geta(𝐶[1])]]
    refnode = getb(𝐶[1])
    popfirst!(𝐶)
    while length(𝐶) > 0
        nextedgeidx = findall(e->geta(e)==refnode, 𝐶)

        if isempty(nextedgeidx)
            push!(regions[numregions], refnode)
            # New region
            numregions += 1
            push!(regions, [geta(𝐶[1])])
            refnode = getb(𝐶[1])
            popfirst!(𝐶)
        else
            if length(nextedgeidx) > 1
                prevnode = regions[numregions][end]
                tempnodes = getb.(𝐶[nextedgeidx])
                idx = findnextnode(prevnode, refnode, tempnodes, geom2fcn)
                nextedgeidx = nextedgeidx[idx]
            else
                nextedgeidx = nextedgeidx[1]
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
"""
function rootsandpoles(
    regions::Vector{Vector{IndexablePoint2D}},
    quadrants::Vector{UInt8},
    geom2fcn::Function
    )

    numregions = size(regions, 1)
    q = Vector{Union{Missing, Int}}(undef, numregions)
    z = Vector{ComplexF64}(undef, numregions)
    for ii in eachindex(regions)
        # XXX: ORDER OF REGIONS??? XXX
        quadrantsequence = [convert(Int8, quadrants[getindex(node)]) for node in regions[ii]]
        # Sign flip because `regions[ii]` are in opposite order of Matlab??
        dQ = -diff(quadrantsequence)
        for jj in eachindex(dQ)
            dQ[jj] == 3 && (dQ[jj] = -1)
            dQ[jj] == -3 && (dQ[jj] = 1)
            # ``|ΔQ| = 2`` is ambiguous; cannot tell whether phase increases or decreases by two quadrants
            abs(dQ[jj]) == 2 && (dQ[jj] = 0)
        end
        q[ii] = sum(dQ)/4
        z[ii] = mean(geom2fcn.(regions[ii]))
    end
    zroots = [z[i] for i in eachindex(z) if q[i] > 0]
    zroots_multiplicity = filter(x->x>0, q)

    zpoles = [z[i] for i in eachindex(z) if q[i] < 0]
    zpoles_multiplicity = filter(x->x<0, q)

    return zroots, zroots_multiplicity, zpoles, zpoles_multiplicity
end

"""
    tesselate!(tess, newnodes, fcn, geom2fcn, tolerance)
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
    quadrants = Vector{UInt8}()

    iteration = 0
    while (iteration < maxiterations) & (numnodes < maxnodes)
        iteration += 1

        # Determine which quadrant function value belongs at each node
        numnewnodes = length(newnodes)
        append!(quadrants, Vector{UInt8}(undef, numnewnodes))
        assignquadrants!(quadrants, newnodes, fcn)

        # Add new nodes to `tess`
        push!(tess, newnodes)
        numnodes += numnewnodes

        # Determine candidate edges that may be near a root or pole
        𝓔, phasediffs = candidateedges(tess, quadrants)
        isempty(𝓔) && error("No roots in the domain")

        # Select candidate edges that are longer than the chosen tolerance
        select𝓔 = filter(e -> longedge(e, tolerance, geom2fcn), 𝓔)
        isempty(select𝓔) && return tess, 𝓔, quadrants, phasediffs
        max𝓔length = maximum(distance(p1, p2) for (p1, p2) in geom2fcn.(select𝓔))
        @debug "Candidate edges length max: $max𝓔length"
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
    grpf(tess, newnodes, fcn, geom2fcn, tolerance)
"""
function grpf(
    tess::DelaunayTessellation2D{IndexablePoint2D},
    newnodes::Vector{IndexablePoint2D},
    fcn::Function,
    geom2fcn::Function,
    tolerance
    )

    tess, 𝓔, quadrants = tesselate!(tess, newnodes, fcn, geom2fcn, tolerance)
    𝐶 = contouredges(tess, 𝓔)
    regions = evaluateregions!(𝐶, geom2fcn)
    zroots, zroots_multiplicity, zpoles, zpoles_multiplicity = rootsandpoles(regions, quadrants, geom2fcn)

    return zroots, zroots_multiplicity, zpoles, zpoles_multiplicity
end

end # module
