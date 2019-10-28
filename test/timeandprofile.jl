using LinearAlgebra
using SpecialFunctions
using VoronoiDelaunay

using GRPF

"""
    mapfunctionval(z, ra, rb, ia, ib)

Linearly map function values within domain from `min_coord` to `max_coord`.
"""
function mapfunctionval(z, ra, rb, ia, ib)
    zr = ra*real(z) + rb
    zi = ia*imag(z) + ib
    return complex(zr, zi)
end

"""
    geom2fcn(pt, ra, rb, ia, ib)

Linearly map geometry values ∈ {`min_coord`, `max_coord`} to domain bounds.

Note: There are floating point errors when converting back and forth.
"""
function geom2fcn(pt::IndexablePoint2D, ra, rb, ia, ib)
    return complex((getx(pt) - rb)/ra, (gety(pt) - ib)/ia)
end
function geom2fcn(edge::VoronoiDelaunay.DelaunayEdge{IndexablePoint2D}, ra, rb, ia, ib)
    return geom2fcn(geta(edge), ra, rb, ia, ib), geom2fcn(getb(edge), ra, rb, ia, ib)
end

function wvgd(z)
      ns = 0.065-4im
      n1 = 1.5835
      nc = 1.0
      d1 = 1.81e-6
      λ₀ = 0.6328e-6
      k₀ = 2π/λ₀
      k₀d1 = k₀*d1
      κ1 = sqrt(n1^2 - z^2)
      γs = sqrt(z^2 - ns^2)
      γc = sqrt(z^2 - nc^2)
      m11 = cos(κ1*k₀d1)
      m12 = im/κ1*sin(κ1*k₀d1)
      m21 = im*κ1*sin(κ1*k₀d1)
      m22 = cos(κ1*k₀d1)
      w = det([1.0    -m11+im*γc*m12
               im*γs  -m21+im*γc*m22])
end

# Analysis parameters
xb = 1.0  # real part begin
xe = 2.5  # real part end
yb = -1.0  # imag part begin
ye = 1.0  # imag part end
r = 0.5  # initial mesh step

const tolerance = 1e-9

origcoords = rectangulardomain(complex(xb, yb), complex(xe, ye), r)

rmin, rmax = minimum(real(origcoords)), maximum(real(origcoords))
imin, imax = minimum(imag(origcoords)), maximum(imag(origcoords))

const ra = (max_coord-min_coord)/(rmax-rmin)
const rb = max_coord - ra*rmax

const ia = (max_coord-min_coord)/(imax-imin)
const ib = max_coord - ia*imax

const origcoordsmapped = mapfunctionval.(origcoords, ra, rb, ia, ib)

fcn(pt) = wvgd(geom2fcn(pt, ra, rb, ia, ib))
geom2fcn(e) = geom2fcn(e, ra, rb, ia, ib)


function gtest()
    newnodes = [IndexablePoint2D(real(coord), imag(coord), idx) for (idx, coord) in enumerate(origcoordsmapped)]
    tess = DelaunayTessellation2D{IndexablePoint2D}(5000)
    GRPF.grpf(tess, newnodes, fcn, geom2fcn, tolerance)
end

function ptest()
    newnodes = [IndexablePoint2D(real(coord), imag(coord), idx) for (idx, coord) in enumerate(origcoordsmapped)]
    tess = DelaunayTessellation2D{IndexablePoint2D}(5000)

    tess, 𝓔, quadrants = GRPF.tesselate!(tess, newnodes, fcn, geom2fcn, tolerance)
    # grpf(tess, newnodes, pt -> wvgd(geom2fcn(pt, ra, rb, ia, ib)),
    #      e -> geom2fcn(e, ra, rb, ia, ib), tolerance)

    return tess, 𝓔, quadrants
end




tess, 𝓔, quadrants = GRPF.tesselate!(tess, newnodes, fcn, geom2fcn, tolerance)



𝐶 = GRPF.contouredges(tess, 𝓔)  # ~2 ms, 1MB alloc
regions = GRPF.evaluateregions!(𝐶, geom2fcn)  # < 1 ms, 10 kB alloc
zroots, zpoles = GRPF.rootsandpoles(regions, quadrants, geom2fcn)

function localrootsandpoles(regions, quadrants, geom2fcn)
    numregions = length(regions)
    # q = Vector{Union{Missing, Int}}(undef, numregions)  # TODO: missing?
    # q = Vector{Int}(undef, numregions)
    # z = Vector{ComplexF64}(undef, numregions)
    zroots = Vector{ComplexF64}()
    zpoles = Vector{ComplexF64}()
    for ii in eachindex(regions)
        # XXX: ORDER OF REGIONS???
        # quadrantsequence = [convert(Int8, quadrants[getindex(node)]) for node in regions[ii]]
        quadrantsequence = [quadrants[getindex(node)] for node in regions[ii]]

        # Sign flip because `regions[ii]` are in opposite order of Matlab??
        dQ = -diff(quadrantsequence)
        for jj in eachindex(dQ)
            dQ[jj] == 3 && (dQ[jj] = -1)
            dQ[jj] == -3 && (dQ[jj] = 1)
            # ``|ΔQ| = 2`` is ambiguous; cannot tell whether phase increases or decreases by two quadrants
            abs(dQ[jj]) == 2 && (dQ[jj] = 0)
        end
        # q[ii] = sum(dQ)/4  # TODO Type?
        # z[ii] = mean(geom2fcn.(regions[ii]))

        q = sum(dQ)/4
        z = mean(geom2fcn.(regions[ii]))
        if q > 0
            push!(zroots, z)
        elseif q < 0
            push!(zpoles, z)
        end
    end

    # for i in eachindex(z)
    #     if q[i] > 0
    #         push!(zroots, z[i])
    #     elseif q[i] < 0
    #         push!(zpoles, z[i])
    #     end
    # end

    return zroots, zpoles
end












numnodes = tess._total_points_added

𝓔 = Vector{GRPF.DelaunayEdge{GRPF.IndexablePoint2D}}()
quadrants = Vector{Int8}()

iteration = 1

# Determine which quadrant function value belongs at each node
numnewnodes = length(newnodes)
append!(quadrants, Vector{Int8}(undef, numnewnodes))
GRPF.assignquadrants!(quadrants, newnodes, fcn)

# Add new nodes to `tess`
push!(tess, newnodes)
numnodes += numnewnodes

# Determine candidate edges that may be near a root or pole
𝓔 = GRPF.candidateedges(tess, quadrants)
isempty(𝓔) && error("No roots in the domain")

# Select candidate edges that are longer than the chosen tolerance
select𝓔 = filter(e -> GRPF.longedge(e, tolerance, geom2fcn), 𝓔)
isempty(select𝓔) && return tess, 𝓔, quadrants

max𝓔length = maximum(GRPF.distance(geom2fcnf(e)) for e in select𝓔)

# How many times does each triangle contain a `select𝓔` node?
trianglecounts = GRPF.counttriangleswithnodes(tess, select𝓔)
zone1triangles, zone2triangles = GRPF.splittriangles(tess, trianglecounts)

# Add new nodes in zone 1
newnodes = Vector{IndexablePoint2D}()
GRPF.zone1newnodes!(newnodes, zone1triangles, geom2fcnf, tolerance)

# Add new nodes in zone 2
GRPF.zone2newnodes!(newnodes, zone2triangles)

# Have to assign indexes to new nodes (which are all currently -1)
setindex!.(newnodes, (1:length(newnodes)).+numnodes)
