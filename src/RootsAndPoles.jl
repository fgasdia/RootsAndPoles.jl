"""
    RootsAndPoles.jl

`RootsAndPoles` is a Julia implementation of the Global complex Roots and Poles Finding
(GRPF) algorithm and the Self-Adaptive Mesh Generator.

Matlab code is available under MIT license at https://github.com/PioKow/GRPF.

# References

[^1]: P. Kowalczyk, “Global complex roots and poles finding algorithm based on phase
    analysis for propagation and radiation problems,” IEEE Transactions on Antennas and
    Propagation, vol. 66, no. 12, pp. 7198–7205, Dec. 2018, doi: 10.1109/TAP.2018.2869213.
[^2]: S. Dziedziewicz, M. Warecka, R. Lech, and P. Kowalczyk, “Self-Adaptive Mesh Generator
    for Global Complex Roots and Poles Finding Algorithm,” IEEE Trans. Microwave Theory
    Techn., vol. 71, no. 7, pp. 2854–2863, July 2023, doi: 10.1109/TMTT.2023.3238014.
"""
module RootsAndPoles

using Base.Threads
using LinearAlgebra, Random
using StaticArrays, ChunkSplitters
using EnumX
using DelaunayTriangulation
const DT = DelaunayTriangulation

const TRIANGLETYPE = NTuple{3,Int}
const EDGETYPE = NTuple{2,Int}

"""
    FinderParams

# Fields

- `tol::Float64 = 1e-9`: maximum edge length of the triangulation before returning.
- `maxiters::Int = 1000`: maximum number of refinement iterations before returning.
- `maxadaptivenodes::Int = 2000`: maximum number of nodes in the triangulation
    built from the adaptive mesh algorithm of SA-GPRF before switching to regular GRPF
    refinement.
- `maxnodes::Int = 20000`: maximum number of nodes in the triangulation before returning.
- `numtasks::Int = 1`: spawn `numtasks` tasks to evaluate the user-provided function
    across the triangulation.
"""
@kwdef struct FinderParams
    tol::Float64 = 1e-9
    maxiters::Int = 1000
    maxadaptivenodes::Int = 2000
    maxnodes::Int = 20000
    numtasks::Int = 1
end

struct ReturnMultiplicity end

abstract type SkinnyMode end

"A triangle is skinny if the ratio of its longest to shortest side is greater than `threshold`."
@kwdef struct MaxMinRatio <: SkinnyMode
    threshold::Float64 = 3
end

"A triangle is skinny if the ratio of its longest side (its \"base\") to its height is greater then `threshold`."
@kwdef struct BaseHeightRatio <: SkinnyMode
    threshold::Float64 = 10
end

include("ComplexMesh.jl")
include("utils.jl")
include("coordinate_domains.jl")

@enumx SpecialEdge::Int8 begin
    Candidate = 1
    Extreme
    Skinny
    Gradient
    Normal = 0
end
getcandidateedges(E::Dict) = (e.first for e in E if e.second == SpecialEdge.Candidate)

struct MeshIterations{M<:ComplexMesh}
    mesh::Vector{M}
    splitedges::Vector{Dict{EDGETYPE,SpecialEdge.T}}  # reason for splitting
    refinement_mode::Vector{Symbol}
end
MeshIterations(m::M) where {M} =
    MeshIterations{M}(M[], Vector{Dict{EDGETYPE,SpecialEdge.T}}(), Vector{Symbol}())
push_mesh!(mi::MeshIterations, m) = push!(mi.mesh, m)
push_edges!(mi::MeshIterations, e) = push!(mi.splitedges, e)
push_mode!(mi::MeshIterations, m) = push!(mi.refinement_mode, m)

function Base.show(io::IO, ::MIME"text/plain", m::MeshIterations)
    println(io, "$(length(m.mesh))-element - MeshIterations")
    print(io, "refinement_mode: ", m.refinement_mode)
end

struct PreviousIteration
    gradients::Dict{TRIANGLETYPE,SVector{2,Float64}}
    splitedges::Set{TRIANGLETYPE}  # includes added midpoints
    triangles::Set{TRIANGLETYPE}
end

struct GradientAnalysis
    gradients::Dict{TRIANGLETYPE,SVector{2,Float64}}
    indicators::Vector{Float64}
    remainingedgelengths::Vector{Float64}
end
Base.empty!(ga::GradientAnalysis) = foreach(f->empty!(getfield(ga, f)), fieldnames(GradientAnalysis))

export FinderParams, MaxMinRatio, BaseHeightRatio
export ComplexMesh, MeshIterations, ReturnMultiplicity, SpecialEdge
export rootsandpoles


"""
    quadrant(z)

Return quadrant on the complex plane of complex number `z`.

| Quadrant |       Phase       |
|:--------:|:-----------------:|
|    1     | 0 ≤ arg f < π/2   |
|    2     | π/2 ≤ arg f < π   |
|    3     | π ≤ arg f < 3π/2  |
|    4     | 3π/2 ≤ arg f < 2π |

Return quadrant 2 if `isnan(z)`.
"""
function quadrant(z)
    r, i = reim(z)
    if r > 0 && i >= 0
        return 1
    elseif r <= 0 && i >= 0
        return 2
    elseif r < 0 && i < 0
        return 3
    elseif r >= 0 && i < 0
        return 4
    elseif isnan(z)
        return 2  # z is possibly a pole
    else
        error("quadrant could not be evaluated for $z")
    end
end

"""
    evalfcn!(f, mesh, startindex; numtasks=Threads.nthreads())

Evaluate function `f` from `startindex` to the last node of `mesh`, updating `mesh.fvals` in place.

Function evaluation is multithreaded into number `numtasks` of tasks.
"""
function evalfcn!(f, mesh, startindex; numtasks=Threads.nthreads())
    indices = startindex:lastindex(mesh)
    @sync for inds in chunks(indices; n=numtasks, split=RoundRobin())
        @spawn for i in inds
            z = complex(coord(mesh, i)...)
            v = f(z)
            fval!(mesh, i, v)
        end
    end
end


"""
    findcandidateedges(mesh)

Return a `Set` of edges within `mesh` for which the change in [`quadrant`](@ref) between
nodes is 2.

Roots or poles are located where the regions described by four different quadrants meet.
Since any triangulation of the four nodes located in the four different quadrants requires
at least one edge of ``|ΔQ| = 2``, then all such edges are candidates for locating a
root or pole.
"""
findcandidateedges(mesh) = Set(minmax(e...) for e in each_solid_edge(mesh) if iscandidateedge(e, mesh))

function findcandidateedges!(edgestosplit::Dict, mesh)
    for e in each_solid_edge(mesh)
        if iscandidateedge(e, mesh)
            edgestosplit[minmax(e...)] = SpecialEdge.Candidate
        end 
    end
end

function iscandidateedge(edge, mesh)
    A, B = edge_vertices(edge)
    qA, qB = quadrant(mesh, A, B)
    ΔQ = abs(qA - qB)
    return ΔQ == 2
end

"""
    findextremeedges(splitedges, mesh)

Return edges for further refinement if the phase of nodes added in the last iteration is
is outside the range of the phases of the original edge nodes.

For node ``m`` added between nodes ``a`` and ``b``, if ``min{ϕₐ,ϕᵦ} ≤ ϕₘ ≤ min{ϕₐ,ϕᵦ}`` is
not true, then every edge between ``m`` and its neighboring nodes should be split.

# References

See Figure 5 of Dziedziewicz, et al., “Self-Adaptive Mesh Generator for Global Complex Roots
and Poles Finding Algorithm,” doi: 10.1109/TMTT.2023.3238014.
"""
function findextremeedges!(edgestosplit::Dict, splitedges, mesh)
    for e in splitedges
        A, B, M = e
        @assert M > A && M > B "node index M is not greater than A and B"
        if unorderedphase(mesh, A, B, M)  # shouldsplit
            for V in get_neighbours(mesh, M)
                DT.is_ghost_vertex(V) || get!(edgestosplit, minmax(M, V), SpecialEdge.Extreme)
            end
        end
    end
end

function unorderedphase(mesh, A, B, M)
    fA, fB, fM = fval(mesh, A, B, M)
    qA, qB, qM = quadrant(fA), quadrant(fB), quadrant(fM)
    pA, pB, pM = angle(fA), angle(fB), angle(fM)

    # `angle` returns a number -π ≤ angle(z) ≤ π, therefore Q3 and Q4 have ϕ < 0
    # to check pA ≤ pM ≤ pB, we add 2π to phases in Q3 if the other node is in Q2.
    # No need to check Q4 because that would be ΔQ = 2, which makes it a candidate edge.
    #==
              0
        Q4    |    Q1
              |
    -π/2--------------π/2 
              |
        Q3    |    Q2
            -π,π
    ==#

    abs(qA - qB) == 2 && return false

    if qA == 2 && qB == 3
        pB += 2π
        qM == 3 && (pM += 2π)
    elseif qA == 3 && qB == 2
        pA += 2π
        qM == 3 && (pM += 2π)
    end
    
    if (pA <= pM <= pB) || (pA >= pM >= pB)
        Δp = abs(pA - pB)
        return Δp >= π  # should be false, not sure if it's possible to be true
    else
        return true  # is extreme node, must split
    end
end

"""
    findskinnyedges(triangles, mesh, skinny_mode)

Return the longest edge of "skinny" triangles as defined by `refinement_mode`. Returns `nothing` if
there are no skinny triangles in `triangles`.

# Modes

- `:baseheightratio`: triangle is skinny if the ratio of the longest edge to the triangle's
height is greater than `skinnyratio`, typically 10
- `:maxminratio`: triangle is skinny if the ratio of the longest edge to the shortest edge
is greater than `skinnyratio`, typically 3
"""
function findskinnyedges!(edgestosplit::Dict, triangles, mesh, skinny_mode::SkinnyMode)
    for T in triangles
        e = _skinnyedge(T, mesh, skinny_mode)
        isnothing(e) || get!(edgestosplit, minmax(e...), SpecialEdge.Skinny)
    end
end

function _skinnyedge(T, mesh, skinny_mode::SkinnyMode)
    A, B, C = DT.triangle_edges(T)

    lA = edge_length(mesh, A)
    lB = edge_length(mesh, B)
    lC = edge_length(mesh, C)

    maxlength, imax = findmax((lA, lB, lC))
    if isskinny(lA, lB, lC, maxlength, skinny_mode)
        # return the longest edge
        imax == 1 && return A
        imax == 2 && return B
        return C  # imax == 3
    else
        return nothing
    end
end

function isskinny(lA, lB, lC, maxlength, m::MaxMinRatio)
    minlength = minimum((lA, lB, lC))
    return maxlength/minlength > m.threshold
end

function isskinny(lA, lB, lC, maxlength, m::BaseHeightRatio)
    s = (lA + lB + lC)/2
    area = sqrt(abs(s*(s - lA)*(s - lB)*(s - lC)))  # abs in case s ≈ lA, lB, or lC and term is negative
    h = 2*area/maxlength  # area = 1/2 × b × h
    return maxlength/h > m.threshold
end

"Return all edges of `mesh` not in `edgestosplit` with a length of at least `atol`."
function remainingedges!(edgestocheck, mesh, edgestosplit, atol)
    K = keys(edgestosplit)
    @assert all(issorted, K) "each edge of `edgestosplit` must be sorted with `minmax`"
    empty!(edgestocheck)
    for e in each_solid_edge(mesh)
        se = minmax(e...)
        if !in(se, K) && (edge_length(mesh, e) > atol)
            push!(edgestocheck, se)
        end
    end
end

"""
    analyzegradients!(edgestosplit, mesh, edges, prev_gradients, minlength, numtosplit)

Return up to number `numtosplit` of `edges` for mesh refinement based on indicator value
``I = α log₁₀(𝓁/minlength)`` where ``α`` is the angle between phase gradient vectors for
triangles on each side of an edge and ``𝓁`` is the edge length. `minlength` is the length of
the shortest edge in the mesh.

See also: [`trianglegradient`](@ref)
"""
function analyzegradients!(edgestosplit::Dict, ga::GradientAnalysis, mesh, edges::Vector,
    prev_gradients, minlength, numtosplit)
    @assert all(issorted, edges) "each edge of `edges` must be sorted by `minmax`"

    # Split every edge if there are fewer edges than we need to split
    numedges = length(edges)
    if numedges <= numtosplit
        foreach(e -> get!(edgestosplit, e, SpecialEdge.Gradient), edges)
        return nothing
    end

    append!(ga.indicators, 
        edgeindicator!(ga.gradients, e, mesh, prev_gradients, minlength) for e in edges)
    p = partialsortperm(ga.indicators, 1:numtosplit, rev=true)
    numposindicators = count(>(eps(eltype(ga.indicators))), ga.indicators)

    if numposindicators < numtosplit
        if numposindicators > 0
            ve = view(edges, view(p, 1:numposindicators))
            foreach(e -> get!(edgestosplit, e, SpecialEdge.Gradient), ve)
            numtosplit -= numposindicators
        end

        # Fill the remaining by longest edge first up to the remaining numtosplit
        zerop = view(p, numposindicators+1:lastindex(p))
        append!(ga.remainingedgelengths, edge_length(mesh, e) for e in view(edges, zerop))
        l = partialsortperm(ga.remainingedgelengths, 1:numtosplit, rev=true)
        ve = view(edges, l)
        foreach(e -> get!(edgestosplit, e, SpecialEdge.Gradient), ve)
    else  # numposindicators >= numtosplit
        ve = view(edges, p)
        foreach(e -> get!(edgestosplit, e, SpecialEdge.Gradient), ve)
    end
    return nothing
end

"Edge \"indicator\" where larger values signify better candidates for splitting."
function edgeindicator!(gradients, e, mesh, prev_gradients, minlength)
    u = get_adjacent(mesh, e)
    v = get_adjacent(mesh, e[2], e[1])  # reverse edge for triangle on other side

    # Edge is on boundary
    if DT.is_ghost_vertex(u) || DT.is_ghost_vertex(v)
        return zero(Float64)
    end

    Tu = (e[1], e[2], u)
    Tv = (e[2], e[1], v)

    gu = gradients[Tu] = get(prev_gradients, Tu, trianglegradient(Tu, mesh))
    gv = gradients[Tv] = get(prev_gradients, Tv, trianglegradient(Tv, mesh))

    # Calculate angle `α` between phase gradients of adjacent triangles (that share an edge)
    α = atan(norm(cross(gu, gv)), dot(gu, gv))
    isnan(α) && (α = π)

    d = edge_length(mesh, e)
    I = α*log10(d/minlength)

    return I
end

"""
    trianglegradient(T, mesh)

Return phase gradient vector for triangle `T`.

# References

See Appendix A of Dziedziewicz, et al., “Self-Adaptive Mesh Generator for Global Complex Roots
and Poles Finding Algorithm,” doi: 10.1109/TMTT.2023.3238014.
"""
function trianglegradient(T, mesh)
    A, B, C = triangle_vertices(T)
    
    pA, pB, pC = angle.(fval(mesh, A, B, C))

    # Phase differences
    ΔpAB = phasediff(pA, pB)
    ΔpBC = phasediff(pB, pC)
    ΔpCA = phasediff(pC, pA)

    circulation = ΔpAB + ΔpBC + ΔpCA

    # Gradient can only be unambiguously determined if the phases sum to zero
    # this is not the case if the triangle has a candidate edge, for example
    if abs(circulation) > 1e-12
        g = SVector(0.0, 0.0)
    else
        cA, cB, cC = coord(mesh, A, B, C)

        W = @SMatrix [cB[1]-cA[1] cB[2]-cA[2]
                      cC[1]-cA[1] cC[2]-cA[2]]
        Wd = SVector(ΔpAB, -ΔpCA)

        g = W\Wd
    end

    return g
end

"Maintain -π ≤ Δp ≤ π"
function phasediff(pA, pB)
    Δp =  pB - pA
    if Δp > π
        Δp -= 2π
    elseif Δp < -π
        Δp += 2π
    end
    return Δp
end

"""
    midpointcoords!(newcoords, splitedges, mesh, edgestosplit)

Update `newcoords` in place with coordinates at the midpoint of each edge in `edgestosplit`.
Also updates `splitedges` in place with tuples of `(u, v, i)` for vertices `u, v` of each
edge and what will be the index of the new vertex in the mesh, `i`.
"""
function midpointcoords!(splitedges, mesh, edgestosplit)
    empty!(splitedges)
    offset = num_solid_vertices(mesh)
    for (i, e) in enumerate(keys(edgestosplit))
        u, v = edge_vertices(e)
        mpt = midpoint(mesh, u, v)

        # Must check if mpt is already in mesh.coords, otherwise a degenerate triangle may
        # cause DelaunayTriangulation to hang
        if !(mpt in mesh.coords)
            push!(mesh.coords, midpoint(mesh, u, v))
            push!(splitedges, (u, v, offset+i))
        end
    end
end

"""
    adaptive()
"""
function adaptive!(edgestosplit::Dict, ga::GradientAnalysis, mesh, previter, edgestocheck,
    tol, skinny_mode::SkinnyMode=BaseHeightRatio())

    # Set edge length tolerance for adaptive refinement
    minlength = minimum(e->edge_length(mesh, e), each_solid_edge(mesh))
    atol = max(minlength, tol)

    # Check condition (5) of SAGRPF for new node zM added to the center of zA and zB:
    # min{ϕA, ϕB} ≤ ϕM ≤ max{ϕA, ϕB}
    # If condition is not met, split all edges attached to node zM
    findextremeedges!(edgestosplit, previter.splitedges, mesh)

    # Always split at least one edge in gradient analysis. Don't count skinny edges.
    numtosplit = max(1, length(edgestosplit))

    # Longest edge of "skinny" triangles
    triangles = (T for T in each_solid_triangle(mesh) if !in(T, previter.triangles))
    findskinnyedges!(edgestosplit, triangles, mesh, skinny_mode)

    remainingedges!(edgestocheck, mesh, edgestosplit, atol)
    filter!(e->edge_length(mesh, e.first) > tol, edgestosplit)
    analyzegradients!(edgestosplit, ga, mesh, edgestocheck, previter.gradients, minlength, numtosplit)

    copy!(previter.gradients, ga.gradients)
    empty!(previter.triangles)  # avoids allocating a Set of each_solid_triangle in copy!
    union!(previter.triangles, each_solid_triangle(mesh))

    return nothing
end

"Return `e` if `e` is in `k` or return `reverse(e)` if `reverse(e)` is in `k`. Otherwise, return `nothing`."
function edgekey(k, e)
    re = DT.reverse_edge(e)
    for v in k
        v == e && return e
        v == re && return re
    end
    return nothing
end

function increment_solid!(d, e, k)
    DT.is_ghost_edge(e) && return nothing
    ee = edgekey(k, e)  # existing edge
    isnothing(ee) ? d[e] = 1 : d[ee] += 1
    return nothing
end

findcandidatetriangles(mesh, candidateedges) =
    Set(Iterators.flatten(candidatetriangle(mesh, e) for e in candidateedges))

function candidatetriangle(mesh, e)
    # Need to gather candidate triangles first, rather than edges
    # Otherwise, two candidate edges part of the same triangle will cause an external
    # edge to be counted twice and deleted when it shouldn't be
    u = get_adjacent(mesh, e)
    v = get_adjacent(mesh, e[2], e[1])  # reverse edge for triangle on other side

    return DT.sort_triangle(e[1], e[2], u), DT.sort_triangle(e[2], e[1], v)
end

function contouredges(triangles)
    edgecount = Dict{EDGETYPE,Int8}()
    for T in triangles
        k = keys(edgecount)
        A, B, C = DT.triangle_edges(T)
        foreach(ee -> increment_solid!(edgecount, ee, k), (A, B, C))
    end
    return edgecount
end

function refinemesh!(f, mesh; params::FinderParams, refinement_mode=:adaptive,
    skinny_mode::SkinnyMode=BaseHeightRatio(),
    iterations::Union{Nothing,MeshIterations}=nothing)

    (; maxiters, maxnodes, maxadaptivenodes, tol, numtasks) = params

    previter = PreviousIteration(
        Dict{TRIANGLETYPE,SVector{2,Float64}}(),  # gradients
        Set{TRIANGLETYPE}(),  # prev_splitedges
        Set{TRIANGLETYPE}(),  # prev_triangles
    )

    ga = GradientAnalysis(
        Dict{TRIANGLETYPE,SVector{2,Float64}}(),
        Vector{Float64}(),
        Vector{Float64}()
    )

    edgestosplit = Dict{EDGETYPE, SpecialEdge.T}()
    edgestocheck = Vector{EDGETYPE}()  # ordered Vector needed
    candidatetriangles = Set{TRIANGLETYPE}()

    iter = 0
    startind = 1
    while iter < maxiters && num_solid_vertices(mesh) < maxnodes && refinement_mode in (:adaptive, :regular)
        iter += 1

        if iter > 1
            startind = num_solid_vertices(mesh) + 1
            addpoints!(mesh)
        end

        evalfcn!(f, mesh, startind; numtasks)
        empty!(edgestosplit)
        findcandidateedges!(edgestosplit, mesh)

        if refinement_mode == :adaptive
            empty!(ga)
            adaptive!(edgestosplit, ga, mesh, previter, edgestocheck, tol, skinny_mode)

            if num_solid_vertices(mesh) >= maxadaptivenodes
                refinement_mode = :regular
            end
        end

        if refinement_mode == :regular
            empty!(candidatetriangles)
            candidateedges = getcandidateedges(edgestosplit)
            union!(candidatetriangles, Iterators.flatten(candidatetriangle(mesh, e) for e in candidateedges))
            findskinnyedges!(edgestosplit, candidatetriangles, mesh, skinny_mode)
            filter!(e -> edge_length(mesh, e.first) > tol, edgestosplit)
        end

        candidateedges = getcandidateedges(edgestosplit)
        if !isempty(candidateedges) && all(e -> edge_length(mesh, e) < tol, candidateedges)
            refinement_mode = :accuracy_achieved
        end

        if isempty(edgestosplit)
            refinement_mode = :aborted_empty_edgestosplit
        else
            midpointcoords!(previter.splitedges, mesh, edgestosplit)
        end

        if !isnothing(iterations)
            push_mesh!(iterations, copy(mesh))
            push_edges!(iterations, copy(edgestosplit))
            push_mode!(iterations, refinement_mode)
        end
    end

    if iter >= maxiters
        @info "Refinement aborted. `maxiters` reached."
    elseif num_solid_vertices(mesh) >= maxnodes
        candidateedges = getcandidateedges(edgestosplit)
        @info "Refinement aborted. `maxnodes` reached."
    end
end

function buildregions(mesh)
    candidateedges = findcandidateedges(mesh)
    triangles = findcandidatetriangles(mesh, candidateedges)
    edgecount = contouredges(triangles)
    filter!(p->p.second == 1, edgecount)  # select edges that were part of a single triangle

    contours = Vector{Set{EDGETYPE}}()
    contourcount = 1

    # Pick an edge from edgecount and then look for a connecting edge among the remaining edges
    k, _ = pop!(edgecount)
    push!(contours, Set((k,)))  # must be done for every new contour
    a, b = k  # vertex `b` is starting point

    while !isempty(edgecount)
        closedcontour = false
        matchfound = false
        for k in keys(edgecount)  # edges
            u, v = k

            if u == b
                delete!(edgecount, k)
                push!(contours[contourcount], k)
                v == a && (closedcontour = true)
                b = v
                matchfound = true
            elseif v == b
                delete!(edgecount, k)
                push!(contours[contourcount], k)
                u == a && (closedcontour = true)
                b = u
                matchfound = true
            end

            if closedcontour && !isempty(edgecount)
                #==
                Edge case: Two contours share a node. `closedcontour` is satisifed when one
                of them closes. This is not desired - they are likely a single root/pole.
                  ________
                 /        \  _____
                /          \/     |
                \          .      |
                 \________/\_____/

                ==#
                k = tangentcontour(edgecount, contours[contourcount])
                if isnothing(k)
                    contourcount += 1
                    k, _ = pop!(edgecount)
                    push!(contours, Set((k,)))
                    a, b = k
                    break
                else
                    delete!(edgecount, k)
                    push!(contours[contourcount], k)
                    a, b = k
                    break
                end
            elseif matchfound
                # keys(edgecount) has changed, so break out of the for loop (and reenter)
                break
            end
        end
        if !matchfound
            contourcount += 1
            k, _ = pop!(edgecount)
            push!(contours, Set((k,)))  # must be done for every new contour
            a, b = k  # vertex `b` is starting point
        end
    end
    return contours
end

"If a vertex of `contour` is still in `edgecount`, return `(a, b)` where `a` is the tangent."
function tangentcontour(edgecount, contour)
    for k in keys(edgecount)
        a, b = k
        for c in contour
            u, v = c
            if a == u || a == v
                return (a, b)
            elseif b == u || b == v
                return (b, a)
            end
        end
    end
    return nothing
end

"Average location of contour nodes."
function centerofmass(mesh, contour)
    p = Set(coord(mesh, c) for c in Iterators.flatten(contour))  # Set removes redundant nodes
    z = sum(x->complex(x...), p)/length(p)  # average of contour nodes
    return z
end

function evaluateregion(mesh, contour)
    q = zero(Float64)  # gets divided by 4
    for e in contour
        A, B = edge_vertices(e)
        
        # Shortcut in case one of the vertices is exactly a root or a pole
        zA, zB = fval(mesh, A, B)
        if iszero(zA)
            return (complex(coord(mesh, A)...), 1)
        elseif iszero(zB)
            return (complex(coord(mesh, B)...), 1)
        elseif isinf(zA) || isnan(zA)
            return (complex(coord(mesh, A)...), -1)
        elseif isinf(zB) || isnan(zB)
            return (complex(coord(mesh, B)...), -1)
        end

        qA, qB = quadrant(mesh, A, B)
        ΔQ = qB - qA
        ΔQ == 3 && (ΔQ = -1)
        ΔQ == -3 && (ΔQ = 1)
        abs(ΔQ) == 2 && @warn "|ΔQ| == 2 should not be along candidate boundary"
        q += ΔQ
    end
    q /= 4
    z = centerofmass(mesh, contour)

    return z, q
end

function evaluateregions(mesh, contours)
    roots = Vector{ComplexF64}()
    poles = similar(roots)
    for c in contours
        z, q = evaluateregion(mesh, c)
        if isroot(q)
            push!(roots, z)
        elseif ispole(q)
            push!(poles, z)
        end
    end
    return roots, poles
end

function evaluateregions(mesh, contours, ::ReturnMultiplicity)
    roots = Vector{ComplexF64}()
    poles = similar(roots)
    rootmultiplicity = Vector{Float64}()
    polemultiplicity = similar(rootmultiplicity)
    for c in contours
        z, q = evaluateregion(mesh, c)
        if isroot(q)
            push!(roots, z)
            push!(rootmultiplicity, q)
        elseif ispole(q)
            push!(poles, z)
            push!(polemultiplicity, q)
        end
    end
    return roots, poles, rootmultiplicity, polemultiplicity
end

isroot(q) = q > 0
ispole(q) = q < 0

"""
    rootsandpoles(f, mesh::ComplexMesh;
        params=FinderParams(),
        refinement_mode=:adaptive,
        dedupe=true,
        skinny_mode::SkinnyMode=BaseHeightRatio(),
        iterations::Union{Nothing,MeshIterations}=nothing)

Return `roots` and `poles` of the function `f` given the initial complex plane `mesh` over
which to search.

- `dedupe` removes duplicate roots and poles within `params.tol`.
- `iterations` saves results at each iteration of the mesh refinement.

!!! warning

    `mesh` is modified in place.
"""
function rootsandpoles(f, mesh::ComplexMesh;
    params=FinderParams(),
    refinement_mode=:adaptive,
    dedupe=true,
    skinny_mode::SkinnyMode=BaseHeightRatio(),
    iterations::Union{Nothing,MeshIterations}=nothing)

    refinemesh!(f, mesh; params, refinement_mode, skinny_mode, iterations)

    contours = buildregions(mesh)
    roots, poles = evaluateregions(mesh, contours)

    if dedupe
        dedupe!(roots, poles; tol=params.tol)
    end

    return roots, poles
end

"""
    rootsandpoles(f, mesh::ComplexMesh, ::ReturnMultiplicity;
        params=FinderParams(),
        refinement_mode=:adaptive,
        dedupe=true,
        skinny_mode::SkinnyMode=BaseHeightRatio(),
        iterations::Union{Nothing,MeshIterations}=nothing)

Return `roots`, `poles`, `rootmult`, and `polemult` multiplicities if `ReturnMultiplicity()`
is also provided to the `rootsandpoles` function.
"""
function rootsandpoles(f, mesh::ComplexMesh, ::ReturnMultiplicity;
    params=FinderParams(),
    refinement_mode=:adaptive,
    dedupe=true,
    skinny_mode::SkinnyMode=BaseHeightRatio(),
    iterations::Union{Nothing,MeshIterations}=nothing)

    refinemesh!(f, mesh; params, refinement_mode, skinny_mode, iterations)

    contours = buildregions(mesh)
    roots, poles, rootmult, polemult = evaluateregions(mesh, contours, ReturnMultiplicity())
    
    if dedupe
        oroots = copy(roots)
        opoles = copy(poles)
        dedupe!(roots, poles; tol=params.tol)

        # Remove multiplicities that had their roots/poles deduped
        Ir = findall(in(roots), oroots)
        Ip = findall(in(poles), opoles)

        keepat!(rootmult, Ir)
        keepat!(polemult, Ip)
    end

    return roots, poles, rootmult, polemult
end

end  # module
