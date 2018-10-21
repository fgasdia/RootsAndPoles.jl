"""
# GRPF: Global complex Roots and Poles Finding algorithm

A Julia implementation of the [Matlab code](https://github.com/PioKow/GRPF) by Piotr Kowalczyk.
"""
# module GRPF

using LinearAlgebra
using Statistics
using StaticArrays
import GeometricalPredicates: intriangle
using VoronoiDelaunay
include("VoronoiDelaunayExtensions.jl")
include("GeneralFunctions.jl")

const maxiterations = 100
const maxnodes = 500000
const skinnytriangle = 3


"""
Generate initial mesh node coordinates for a rectangular domain ∈ {`Zb` to `Ze`}.

See also: `rect_dom.m`
"""
function rectangulardomain(Zb, Ze, Δr)
    X = real(Ze) - real(Zb)
    Y = imag(Ze) - imag(Zb)

    n = ceil(Int64, Y/Δr + 1)
    dy = Y/(n-1)
    m = ceil(Int64, X/sqrt(Δr^2 - dy^2/4) + 1)
    dx = X/(m-1)

    vx = range(real(Zb), length=m, stop=real(Ze))
    vy = range(imag(Zb), length=n, stop=imag(Ze))

    tmp = ones(n)
    tmp[n] = 0.0

    y = repeat(vy, m)
    y .+= 0.5*dy*kron((1 .+ (-1).^(1:m))/2, tmp)
    y = [y; fill(imag(Zb), size(2:2:m))]
    x = reshape(repeat(vx', n), m*n)
    x = [x; ((2:2:m) .- 1)*dx .+ real(Zb)]

    return complex.(x, y)
end

"""
Generate initial mesh coordinates for a circular disk domain of radius `R` and center (0, 0)
for ``|z| < R``.

See also: `disk_dom.m`
"""
function diskdomain(R, r)
    h = r*sqrt(3)/2
    n = 1 + round(Int64, R/h)
    Rn = (1:n)*R/n
    newnodes = [complex(0.0)]
    f₀ = 0.
    np = 6
    for ii = 1:n
        f = f₀ .+ range(0, stop=2π, length=np+1)
        xyn = Rn[ii]*complex.(cos.(f[1:end-1]), sin.(f[1:end-1]))
        append!(newnodes, xyn)
        f₀ += π/6/n
        np += 6
    end
    return newnodes
end

"""
Converts complex function value `val` to quadrant number.

See also: `vinq.m`
"""
function quadrant(val::ComplexF64)
    if (real(val) > 0) & (imag(val) >= 0)
        quad = 1
    elseif (real(val) <= 0) & (imag(val) > 0)
        quad = 2
    elseif (real(val) < 0) & (imag(val) <= 0)
        quad = 3
    elseif (real(val) >= 0) & (imag(val) < 0)
        quad = 4
    else
        error("Function value cannot be assigned to quadrant.")
    end
end

"""
Evaluate function `fcn` for [`quadrant`](@ref) at `newnodescoords` and fill `quadrants`.

`quadrants` is a Vector{} where each index corresponds to `node` index.
"""
function assignquadrants!(quadrants, nodes, fcn)
    for ii in eachindex(nodes)
        val = fcn(nodes[ii])
        quadrants[getindex(nodes[ii])] = quadrant(val)
    end
end

"""
Return candidate edges `𝓔`.

Since any triangulation of the four nodes located in the four different quadrants requires
at least one edge of ``|ΔQ| = 2``, then all such edges are potentially in the vicinity of a
root or pole.

`phasediffs` is returned only for diagnosis and plotting. It can be removed without affecting
this function.

Note: Order of `𝓔` is not guaranteed.
"""
function candidateedges(tess, quadrants)
    phasediffs = Vector{Int64}()
    𝓔 = Vector{DelaunayEdge}()
    for edge in delaunayedges(tess)
        nodea, nodeb = geta(edge), getb(edge)
        idxa, idxb = getindex(nodea), getindex(nodeb)

        # To match Matlab, force `idxa` < `idxb`
        # (order doesn't matter for `ΔQ == 2`, which is the only case we care about)
        # if idxa > idxb
        #     idxa, idxb = idxb, idxa
        # end

        ΔQ = mod(quadrants[idxa] - quadrants[idxb], 4)
        ΔQ == 2 && push!(𝓔, edge)

        push!(phasediffs, ΔQ)
    end
    return 𝓔, phasediffs
end

"""
Select edges with length greater than `tolerance`.
"""
function selectedges(edges, tolerance, geom2fcn)
    min𝓔length = typemax(Float64)
    max𝓔length = typemin(Float64)
    select𝓔 = Vector{DelaunayEdge}()
    for ii in eachindex(edges)
        elen = distance(geom2fcn(edges[ii])...)[1]  # unpacks single-element array
        elen > tolerance && push!(select𝓔, edges[ii])
        if elen > max𝓔length
            max𝓔length = elen
        elseif elen < min𝓔length
            min𝓔length = elen
        end
    end
    return select𝓔, min𝓔length, max𝓔length
end

"""
Counts how many times each triangle contains a node.
"""
function counttriangleswithnodes(tess, edges::Array{DelaunayEdge,1})
    # Nodes of select edges
    edgenodes = Vector{IndexablePoint2D}()
    uniquenodes!(edgenodes, edges)

    # TODO: This part of function is incredibly slow!
    trianglecounts = zeros(Int64, count(.!isexternal.(tess._trigs)))

    triidx = 0
    for triangle in tess
        triidx += 1
        # `triangle` is in general not equal to `tess._trigs[triidx]`
        ea = geta(triangle)
        eb = getb(triangle)
        ec = getc(triangle)
        for nodeidx in eachindex(edgenodes)
            if (ea == edgenodes[nodeidx]) | (eb == edgenodes[nodeidx]) | (ec == edgenodes[nodeidx])
                trianglecounts[triidx] += 1
            end
        end
    end
    return trianglecounts
end
function uniquenodes!(edgenodes, edges::Array{DelaunayEdge,1})
    for ii in eachindex(edges)
        nodea = geta(edges[ii])
        nodeb = getb(edges[ii])

        nodea in edgenodes || push!(edgenodes, nodea)
        nodeb in edgenodes || push!(edgenodes, nodeb)
    end
    nothing
end


"""
Add nodes to `newnodes` in zone 1.
"""
function zone1newnodes!(newnodes::AbstractArray{IndexablePoint2D,1},
                        triangles::AbstractArray{DelaunayTriangle,1}, geom2fcn, tolerance)
    n1a = geta(triangles[1])
    n1b = getb(triangles[1])
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
function addnewnode!(newnodes::AbstractArray{IndexablePoint2D}, node1::IndexablePoint2D,
                     node2::IndexablePoint2D, geom2fcn, tolerance)
    if distance(geom2fcn(node1), geom2fcn(node2)) > tolerance
        avgnode = (node1+node2)/2
        for ii in eachindex(newnodes)
            distance(newnodes[ii], avgnode) < 2eps() && return nothing
        end
        push!(newnodes, avgnode)  # only executed if we haven't returned
    end
    nothing
end


"""
Add nodes to `newnodes` in zone 2 (skinny triangles).
"""
function zone2newnodes!(newnodes, triangles)
    for ii in eachindex(triangles)
        na = geta(triangles[ii])
        nb = getb(triangles[ii])
        nc = getc(triangles[ii])
        l1 = distance(nb, na)
        l2 = distance(nc, nb)
        l3 = distance(na, nc)
        if max(l1,l2,l3)/min(l1,l2,l3) > skinnytriangle
            avgnode = (na+nb+nc)/3
            push!(newnodes, avgnode)
        end
    end
    nothing
end


"""
Find contour edges from all candidate edges.
"""
function contouredges(tess, edges)
    # Edges of triangles that contain at least 1 of `edges`
    tmpedges = Vector{DelaunayEdge}()
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
    𝐶 = Vector{DelaunayEdge}()
    duplicateedges = zeros(Int64, size(tmpedges, 1))
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
"""
function splittriangles(tess, trianglecounts)
    zone1triangles = Vector{DelaunayTriangle}()
    zone2triangles = Vector{DelaunayTriangle}()
    for (ii, triangle) in enumerate(tess)
        if trianglecounts[ii] > 1
            push!(zone1triangles, triangle)
        elseif trianglecounts[ii] == 1
            push!(zone2triangles, triangle)
        end
    end
    return zone1triangles, zone2triangles
end

"""
Find the index of the next node in the candidate region boundary process. The next one (after
the reference) is picked from the fixed set of nodes.

See also: `FindNextNode.m`
"""
function findnextnode(prevnode, refnode, tempnodes, geom2fcn)
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
TODO: Go in reverse (from `end` rather than `1` so we don't popfirst!) or use careful indexing
and don't pop at all

# Note: The nodes of each region are in reverse order compared to Matlab wrt their quadrants?
"""
function evaluateregions!(𝐶, geom2fcn)
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
"""
function rootsandpoles(regions, quadrants, geom2fcn)
    numregions = size(regions, 1)
    q = Vector{Int64}(undef, numregions)
    z = Vector{ComplexF64}(undef, numregions)
    for ii in eachindex(regions)
        quadrantsequence = [quadrants[getindex(node)] for node in regions[ii]]
        # Sign flip because `regions[ii]` are in opposite order of Matlab
        dQ = -diff(quadrantsequence)
        for jj in eachindex(dQ)
            dQ[jj] == 3 && (dQ[jj] = -1)
            dQ[jj] == -3 && (dQ[jj] = 1)
            # ``|ΔQ| = 2`` is ambiguous; cannot tell whether phase increases or decreases by two quadrants
            abs(dQ[jj]) == 2 && (dQ[jj] = NaN)
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
"""
function tesselate!(tess, newnodes, fcn, geom2fcn, tolerance)
    # Initialize
    numnodes = tess._total_points_added
    @assert numnodes == 0

    𝓔 = Vector{DelaunayEdge}()
    quadrants = Vector{Int64}()

    iteration = 0
    while (iteration < maxiterations) & (numnodes < maxnodes)
        iteration += 1

        # Determine which quadrant function value belongs at each node
        numnewnodes = length(newnodes)
        append!(quadrants, Vector{Int64}(undef, numnewnodes))
        assignquadrants!(quadrants, newnodes, fcn)

        # Add new nodes to `tess`
        push!(tess, newnodes)
        numnodes += numnewnodes

        # Determine candidate edges that may be near a root or pole
        𝓔, phasediffs = candidateedges(tess, quadrants)
        isempty(𝓔) && error("No roots in the domain")

        # Select candidate edges that are longer than the chosen tolerance
        select𝓔, min𝓔length, max𝓔length = selectedges(𝓔, tolerance, geom2fcn)
        @debug "Candidate edges length min: $min𝓔length, max: $max𝓔length"
        max𝓔length < tolerance && return tess, 𝓔, quadrants

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


"""
"""
function main(tess, newnodes, fcn, geom2fcn, tolerance)
    tess, 𝓔, quadrants = tesselate!(tess, newnodes, fcn, geom2fcn, tolerance)
    𝐶 = contouredges(tess, 𝓔)
    regions = evaluateregions!(𝐶, geom2fcn)
    zroots, zroots_multiplicity, zpoles, zpoles_multiplicity = rootsandpoles(regions, quadrants, e -> geom2fcn(e, ra, rb, ia, ib))

    return zroots, zroots_multiplicity, zpoles, zpoles_multiplicity
end


# end # module
