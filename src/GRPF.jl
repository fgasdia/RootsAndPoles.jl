__precompile__(true)

"""
# GRPF: Global complex Roots and Poles Finding algorithm

A Julia implementation of the GRPF algorithm. Matlab code is available under MIT
license at https://github.com/PioKow/GRPF.

# References

[^1] P. Kowalczyk, “Global complex roots and poles finding algorithm based on phase
analysis for propagation and radiation problems,” IEEE Transactions on Antennas
and Propagation, vol. 66, no. 12, pp. 7198–7205, Dec. 2018,
doi: 10.1109/TAP.2018.2869213.
"""
module GRPF

#==
NOTE: Some variable conversions from the original GRPF papers to this code:

| Paper | Code |
|-------|------|
|   𝓔   |   E  |
|   𝐶   |   C  |
|   ϕ   |  phi |
==#

import Base
using LinearAlgebra
using VoronoiDelaunay
import VoronoiDelaunay: getx, gety, DelaunayEdge, DelaunayTriangle

# NOTE: `max_coord` and `min_coord` are provided by `VoronoiDelaunay.jl`
# We are even more conservative than going from `max_coord` to
# `min_coords` because it is possible to run into floating point issues at
# the very limits
const MAXCOORD = nextfloat(max_coord, -10)
const MINCOORD = nextfloat(min_coord, 10)

"""
    GRPFParams

Structure for holding values used by `GRPF` to stop iterating or split Delaunay
triangles.

`maxiterations` is the maximum number of refinement iterations before `grpf`
returns. By default, `maxiterations` is 100.

`maxnodes` is the maximum number of Delaunay tessalation nodes before `grpf`
returns. By default, `maxnodes` is 500000.

Delaunay triangles with ratio of their longest to shortest side length greater
than `skinnytriangle` will be split during the `grpf` refinement iterations.

`tess_sizehint` is used to provide a size hint to the total number of expected
nodes in the Delaunay tesselation. Setting this number approximately correct
will improve performance. By default, `tess_sizehint` is 5000.

`tolerance` maximum allowed edge length of the tesselation defined in the
`origcoords` domain. By default, `tolerance` is 1e-9.
"""
struct GRPFParams{T<:Real}
    maxiterations::Int
    maxnodes::Int
    skinnytriangle::Int
    tess_sizehint::Int
    tolerance::T
end

"""
    GRPFParams(tess_sizehint::Integer, tolerance::Real)

Convenience function for creating a `GRPFParams` object with the most
important parameters, `tess_sizehint` and `tolerance`.
"""
GRPFParams(tess_sizehint::Integer, tolerance::Real) = GRPFParams(100, 500000, 3, tess_sizehint, tolerance)
GRPFParams() = GRPFParams(100, 500000, 3, 5000, 1e-9)

struct PlotData end

"""
    Geometry2Function

Store conversion coefficients from the `VoronoiDelaunay` domain to the `origcoords`
domain.
"""
struct Geometry2Function{T<:AbstractFloat}
    ra::T
    rb::T
    ia::T
    ib::T
end
(f::Geometry2Function)(z) = geom2fcn(z, f.ra, f.rb, f.ia, f.ib)
(f::Geometry2Function)(x, y) = geom2fcn(x, y, f.ra, f.rb, f.ia, f.ib)
Base.eltype(f::Geometry2Function{T}) where T = T

"""
    ScaledFunction{T<:Function}

Store conversion coefficients from the `VoronoiDelaunay` domain to the
`origcoords` domain so that a `ScaledFunction` can be called in place of the original
function when providing an argument in the `VoronoiDelaunay` domain.
"""
struct ScaledFunction{F<:Function,G<:Geometry2Function}
    fcn::F
    g2f::G
end
(f::ScaledFunction)(z) = f.fcn(f.g2f(z))

# These files need the above structs defined
include("VoronoiDelaunayExtensions.jl")
include("utils.jl")
include("coordinate_domains.jl")

export rectangulardomain, diskdomain, grpf, PlotData, GRPFParams

"""
    quadrant(val)::Int8

Convert complex function value `val` to quadrant number.
"""
@inline function quadrant(val)::Int8
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
    assignquadrants!(quadrants, nodes, f)

Evaluate function `f` for [`quadrant`](@ref) at `nodes` and fill `quadrants`.

`quadrants` is a vector where each index corresponds to `IndexablePoint2D` node
index.
"""
@inline function assignquadrants!(
    quadrants::Vector{<:Integer},
    nodes::Vector{IndexablePoint2D},
    f::ScaledFunction{T}) where T

    for ii in eachindex(nodes)
        p = @inbounds nodes[ii]
        quadrants[getindex(p)] = quadrant(f(p))
    end
    return nothing
end

"""
    candidateedges!(E, tess, quadrants)

Fill in candidate edges `E` that contain a phase change of 2 quadrants.

Any root or pole is located at the point where the regions described by four
different quadrants meet. Since any triangulation of the four nodes located in
the four different quadrants requires at least one edge of ``|ΔQ| = 2``, then
all such edges are potentially in the vicinity of a root or pole.

Order of `E` is not guaranteed.
"""
function candidateedges!(
    E::Vector{DelaunayEdge{IndexablePoint2D}},
    tess::DelaunayTessellation2D{IndexablePoint2D},
    quadrants::Vector{<:Integer}
    )

    # Better performance compared to `delaunayedges()`
    # see: https://github.com/JuliaGeometry/VoronoiDelaunay.jl/issues/47
    @inbounds for ix in 2:tess._last_trig_index
        tr = tess._trigs[ix]
        isexternal(tr) && continue

        # precalculate
        atr, btr, ctr = geta(tr), getb(tr), getc(tr)
        atri, btri, ctri = getindex(atr), getindex(btr), getindex(ctr)

        ix_na = tr._neighbour_a
        if (ix_na > ix) || isexternal(tess._trigs[ix_na])
            ΔQ = mod(quadrants[btri] - quadrants[ctri], 4)  # phase difference
            if ΔQ == 2
                edge = DelaunayEdge(btr, ctr)
                push!(E, edge)
            end
        end

        ix_nb = tr._neighbour_b
        if (ix_nb > ix) || isexternal(tess._trigs[ix_nb])
            ΔQ = mod(quadrants[atri] - quadrants[ctri], 4)
            if ΔQ == 2
                edge = DelaunayEdge(atr, ctr)
                push!(E, edge)
            end
        end

        ix_nc = tr._neighbour_c
        if (ix_nc > ix) || isexternal(tess._trigs[ix_nc])
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
    candidateedges!(E, tess, quadrants, ::PlotData)

Return candidate edges `E` and the `phasediffs` across each edge.
"""
function candidateedges!(
    E::Vector{DelaunayEdge{IndexablePoint2D}},
    tess::DelaunayTessellation2D{IndexablePoint2D},
    quadrants::Vector{<:Integer},
    ::PlotData
    )

    phasediffs = Vector{Int8}()

    for edge in delaunayedges_fast(tess)
        nodea, nodeb = geta(edge), getb(edge)
        idxa, idxb = getindex(nodea), getindex(nodeb)

        # NOTE: To match Matlab, force `idxa` < `idxb`
        # (order doesn't matter for `ΔQ == 2`, which is the only case we care about)
        if idxa > idxb
            idxa, idxb = idxb, idxa
        end

        @inbounds ΔQ = mod(quadrants[idxa] - quadrants[idxb], Int8(4))  # phase difference
        if ΔQ == 2
            push!(E, edge)
        end

        push!(phasediffs, ΔQ)
    end

    return E, phasediffs
end

"""
    zone(triangle, edge_idxs)

Return zone `1` or `2` for `triangle`.
"""
function zone(
    triangle::DelaunayTriangle{IndexablePoint2D},
    edge_idxs::Vector{<:Integer}
    )

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
    uniqueindices(edges)

Return unique indices of all nodes in `edges`.
"""
@inline function uniqueindices(edges::Vector{DelaunayEdge{IndexablePoint2D}})
    idxs = Vector{Int}()
    for i in eachindex(edges)
        @inbounds ei = edges[i]
        push!(idxs, getindex(geta(ei)))
        push!(idxs, getindex(getb(ei)))
    end
    sort!(idxs)  # calling `sort!` first makes `unique!` more efficient
    unique!(idxs)
    return idxs
end

"""
    zone1newnodes!(newnodes, triangles, g2f, tolerance)

Add nodes to `newnodes` in zone 1, i.e. triangles that had more than one node.
"""
function zone1newnodes!(
    newnodes::Vector{IndexablePoint2D},
    triangles::Vector{DelaunayTriangle{IndexablePoint2D}},
    g2f::Geometry2Function,
    tolerance
    )

    triangle1 = triangles[1]
    n1a = geta(triangle1)
    n1b = getb(triangle1)
    push!(newnodes, (n1a+n1b)/2)

    for ii = 1:length(triangles)-1
        @inbounds triangle = triangles[ii]
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

@inline function addnewnode!(
    newnodes::Vector{IndexablePoint2D},
    node1::IndexablePoint2D,
    node2::IndexablePoint2D,
    g2f::Geometry2Function,
    tolerance
    )

    if distance(g2f(node1), g2f(node2)) > tolerance
        avgnode = (node1+node2)/2
        @inbounds for ii in eachindex(newnodes)
            distance(newnodes[ii], avgnode) < 2*eps() && return nothing
        end
        push!(newnodes, avgnode)  # only executed if we haven't already returned
    end
    return nothing
end

"""
    zone2newnodes!(newnodes, triangle, skinnytriangle)

Add node to `newnodes` for zone 2 ("skinny") triangles.

`skinnytriangle` is the maximum allowed ratio of the longest to shortest side
length of each triangle.
"""
@inline function zone2newnode!(
    newnodes::Vector{IndexablePoint2D},
    triangle::DelaunayTriangle{IndexablePoint2D},
    skinnytriangle
    )

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
    splittriangles!(newnodes, tess, edge_idxs, params)

Add zone 2 triangles to `newnodes` and then return a vector of zone 1 triangles,
which require special handling.
"""
function splittriangles!(
    newnodes::AbstractVector,
    tess::DelaunayTessellation2D{IndexablePoint2D},
    edge_idxs::Vector{<:Integer},
    params::GRPFParams
    )

    zone1triangles = Vector{DelaunayTriangle{IndexablePoint2D}}()
    for triangle in tess
        z = zone(triangle, edge_idxs)

        if z == 1
            push!(zone1triangles, triangle)
        elseif z == 2
            zone2newnode!(newnodes, triangle, params.skinnytriangle)
        end
    end
    return zone1triangles
end

"""
    findnextnode(prevnode, refnode, tempnodes, g2f)

Find the index of the next node in the candidate region boundary process. The
next one (after the reference) is picked from the fixed set of nodes.
"""
@inline function findnextnode(
    prevnode::IndexablePoint2D,
    refnode::IndexablePoint2D,
    tempnodes::Vector{IndexablePoint2D},
    g2f::Geometry2Function
    )

    P = g2f(prevnode)
    S = g2f(refnode)

    minphi = 2π + 1  # max diff of angles is 2π, so this is guaranteed larger
    minphi_idx = 1

    for i in eachindex(tempnodes)
        @inbounds N = g2f(tempnodes[i])

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

Find contour edges from all candidate edges.
"""
function contouredges(
    tess::DelaunayTessellation2D{IndexablePoint2D},
    edges::Vector{DelaunayEdge{IndexablePoint2D}}
    )

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
function evaluateregions!(
    C::Vector{DelaunayEdge{IndexablePoint2D}},
    g2f::Geometry2Function
    )

    # NOTE: The nodes of each region are in reverse order compared to Matlab
    # with respect to their quadrants

    # Initialize
    numregions = 1

    regions = [[geta(C[1])]]

    refnode = getb(C[1])  # type annotated to assist with boxing
    popfirst!(C)

    nextedgeidxs = Vector{Int}()
    while length(C) > 0

        # This loop is equivalent to `findall(e->geta(e)==refnode, C)`
        # but avoids closure Core.Box issue
        @inbounds for i in eachindex(C)
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

Identify roots and poles of function based on regions and quadrants.
"""
function rootsandpoles(
    regions::Vector{Vector{IndexablePoint2D}},
    quadrants::Vector{Int8},
    g2f::Geometry2Function{T}
    ) where T

    numregions = length(regions)

    complexT = complex(T)
    zroots = Vector{complexT}()
    zpoles = Vector{complexT}()
    for ii in eachindex(regions)
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
        z = sum(g2f.(regions[ii]))/length(regions[ii])

        if q > 0
            push!(zroots, convert(complexT, z))  # convert in case T isn't Float64
        elseif q < 0
            push!(zpoles, convert(complexT, z))
        end
    end

    return zroots, zpoles
end

"""
    tesselate!(tess, newnodes, fcn, params)

Label quadrants, identify candidate edges, and iteratively split triangles.
"""
function tesselate!(
    tess::DelaunayTessellation2D{IndexablePoint2D},
    newnodes::Vector{IndexablePoint2D},
    f::ScaledFunction{F,G},
    params::GRPFParams
    ) where {F,G}

    # Initialize
    numnodes = tess._total_points_added
    @assert numnodes == 0

    g2f = f.g2f

    E = Vector{DelaunayEdge{IndexablePoint2D}}()
    quadrants = Vector{Int8}()

    iteration = 0
    while (iteration < params.maxiterations) && (numnodes < params.maxnodes)
        iteration += 1

        # Determine which quadrant function value belongs at each node
        numnewnodes = length(newnodes)
        append!(quadrants, Vector{Int8}(undef, numnewnodes))
        assignquadrants!(quadrants, newnodes, f)

        # Add new nodes to `tess`
        push!(tess, newnodes)
        numnodes += numnewnodes

        # Determine candidate edges that may be near a root or pole
        empty!(E)  # start with a blank E
        candidateedges!(E, tess, quadrants)
        isempty(E) && tess, E, quadrants  # no roots or poles found

        # Select candidate edges that are longer than the chosen tolerance
        selectE = filter(e -> longedge(e, params.tolerance, g2f), E)
        isempty(selectE) && return tess, E, quadrants

        # return if maximum edge length has reached tolerance
        maxElength = maximum(distance(g2f(e)) for e in selectE)
        maxElength < params.tolerance && return tess, E, quadrants

        # Get unique indices of nodes in `edges`
        edge_idxs = uniqueindices(selectE)

        # Refine (split) triangles
        newnodes = Vector{IndexablePoint2D}()
        zone1triangles = splittriangles!(newnodes, tess, edge_idxs, params)

        # Add new nodes in zone 1
        zone1newnodes!(newnodes, zone1triangles, g2f, params.tolerance)

        # Have to assign indexes to new nodes (which are all currently -1)
        setindex!.(newnodes, (1:length(newnodes)).+numnodes)
    end

    return tess, E, quadrants
end

function tesselate!(
    tess::DelaunayTessellation2D{IndexablePoint2D},
    newnodes::Vector{IndexablePoint2D},
    f::ScaledFunction{F,G},
    params::GRPFParams,
    ::PlotData
    ) where {F,G}

    # Initialize
    numnodes = tess._total_points_added
    @assert numnodes == 0

    g2f = f.g2f

    E = Vector{DelaunayEdge{IndexablePoint2D}}()
    quadrants = Vector{Int8}()

    iteration = 0
    while (iteration < params.maxiterations) && (numnodes < params.maxnodes)
        iteration += 1

        # Determine which quadrant function value belongs at each node
        numnewnodes = length(newnodes)
        append!(quadrants, Vector{Int8}(undef, numnewnodes))
        assignquadrants!(quadrants, newnodes, f)

        # Add new nodes to `tess`
        push!(tess, newnodes)
        numnodes += numnewnodes

        # Determine candidate edges that may be near a root or pole
        empty!(E)  # always start with a blank E
        E, phasediffs = candidateedges!(E, tess, quadrants, PlotData())
        isempty(E) && return tess, E, quadrants, phasediffs  # no roots or poles found

        # Select candidate edges that are longer than the chosen tolerance
        selectE = filter(e -> longedge(e, params.tolerance, g2f), E)
        isempty(selectE) && return tess, E, quadrants, phasediffs

        # return if maximum edge length has reached tolerance
        maxElength = maximum(distance(g2f(e)) for e in selectE)
        maxElength < params.tolerance && return tess, E, quadrants, phasediffs

        # Get unique indices of nodes in `edges`
        edge_idxs = uniqueindices(selectE)

        # Refine (split) triangles
        newnodes = Vector{IndexablePoint2D}()
        zone1triangles = splittriangles!(newnodes, tess, edge_idxs, params)

        # Add new nodes in zone 1
        zone1newnodes!(newnodes, zone1triangles, g2f, params.tolerance)

        # Have to assign indexes to new nodes (which are all currently -1)
        setindex!.(newnodes, (1:length(newnodes)).+numnodes)
    end

    return tess, E, quadrants, phasediffs
end

"""
    grpf(fcn, origcoords, params=GRPFParams())

Return a vector `roots` and a vector of `poles` of a single (complex) argument
function `fcn`.

Searches within a domain specified by the vector of `origcoords`.

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
function grpf(fcn::Function, origcoords::AbstractArray{<:Complex}, params=GRPFParams())

    # TODO: See pull request #50 on VoronoiDelaunay.jl which handles the space
    # mapping automatically.

    # Need to map space domain for VoronoiDelaunay.jl
    rmin, rmax = minimum(real, origcoords), maximum(real, origcoords)
    imin, imax = minimum(imag, origcoords), maximum(imag, origcoords)

    # Be slightly conservative with our scaling to ensure we stay inside of
    # VoronoiDelaunay.jl `max_coord` and `min_coord`
    width = MAXCOORD - MINCOORD
    ra = width/(rmax-rmin)
    rb = MAXCOORD - ra*rmax

    ia = width/(imax-imin)
    ib = MAXCOORD - ia*imax

    origcoords = fcn2geom.(origcoords, ra, rb, ia, ib)

    @assert minimum(real, origcoords) >= min_coord && minimum(imag, origcoords) >= min_coord &&
        maximum(real, origcoords) <= max_coord && maximum(imag, origcoords) <= max_coord "Scaled coordinates out of bounds"

    newnodes = [IndexablePoint2D(real(coord), imag(coord), idx) for (idx, coord) in enumerate(origcoords)]
    tess = DelaunayTessellation2D{IndexablePoint2D}(params.tess_sizehint)

    g2f = Geometry2Function(ra, rb, ia, ib)
    f = ScaledFunction(fcn, g2f)

    tess, E, quadrants = tesselate!(tess, newnodes, f, params)

    complexT = complex(eltype(g2f))
    isempty(E) && return Vector{complexT}(), Vector{complexT}()

    C = contouredges(tess, E)
    regions = evaluateregions!(C, g2f)
    zroots, zpoles = rootsandpoles(regions, quadrants, g2f)

    return zroots, zpoles
end

"""
    grpf(fcn, origcoords, ::PlotData, params=GRPFParams())

Variant of `grpf` that returns `quadrants` and `phasediffs` in addition to
`zroots` and `zpoles`, primarily for plotting or diagnostics.

# Examples
```jldoctest
julia> simplefcn(z) = (z - 1)*(z - im)^2*(z + 1)^3/(z + im)

julia> xb, xe = -2, 2

julia> yb, ye = -2, 2

julia> r = 0.1

julia> origcoords = rectangulardomain(complex(xb, yb), complex(xe, ye), r)

julia> roots, poles, quadrants, phasediffs = grpf(simplefcn, origcoords, PlotData());
```
"""
function grpf(fcn::Function, origcoords::AbstractArray{<:Complex}, ::PlotData, params=GRPFParams())
    # Need to map space domain for VoronoiDelaunay.jl
    rmin, rmax = minimum(real, origcoords), maximum(real, origcoords)
    imin, imax = minimum(imag, origcoords), maximum(imag, origcoords)

    # Be slightly conservative with our scaling to ensure we stay inside of
    # VoronoiDelaunay.jl `max_coord` and `min_coord`
    width = MAXCOORD - MINCOORD
    ra = width/(rmax-rmin)
    rb = MAXCOORD - ra*rmax

    ia = width/(imax-imin)
    ib = MAXCOORD - ia*imax

    origcoords = fcn2geom.(origcoords, ra, rb, ia, ib)
    @assert minimum(real, origcoords) >= min_coord && minimum(imag, origcoords) >= min_coord &&
        maximum(real, origcoords) <= max_coord && maximum(imag, origcoords) <= max_coord "Scaled coordinates out of bounds"

    newnodes = [IndexablePoint2D(real(coord), imag(coord), idx) for (idx, coord) in enumerate(origcoords)]
    tess = DelaunayTessellation2D{IndexablePoint2D}(params.tess_sizehint)

    g2f = Geometry2Function(ra, rb, ia, ib)
    f = ScaledFunction(fcn, g2f)

    tess, E, quadrants, phasediffs = tesselate!(tess, newnodes, f, params, PlotData())

    complexT = complex(eltype(g2f))
    isempty(E) && return Vector{complexT}(), Vector{complexT}(), quadrants, phasediffs, tess, g2f

    C = contouredges(tess, E)
    regions = evaluateregions!(C, g2f)
    zroots, zpoles = rootsandpoles(regions, quadrants, g2f)

    return zroots, zpoles, quadrants, phasediffs, tess, g2f
end

end # module
