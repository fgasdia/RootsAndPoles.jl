distance(p1::Complex, p2::Complex) = hypot(p2 - p1)
distance((p1, p2)::Tuple{Complex, Complex}) = distance(p1, p2)
function distance(p1::AbstractPoint2D, p2::AbstractPoint2D)
    z1 = complex(p1._x, p1._y)
    z2 = complex(p2._x, p2._y)
    return distance(z1, z2)
end

"""
    longedge(edge, tolerance, geom2fcn)

Return true if `edge` has length greater than `tolerance`.
"""
function longedge(edge::DelaunayEdge, tolerance, g2f::Geometry2Function)
    return distance(g2f(edge)) > tolerance
end

"""
    fcn2geom(x, xmin, xmax)

Linearly scale `x`, which exists ∈ (`xmin`, `xmax`) into `RootsAndPoles` domain.
"""
function fcn2geom(x, xmin, xmax)
    return (MAXCOORD - MINCOORD)*(x - xmin)/(xmax - xmin) + MINCOORD
end

# TODO: See pull request #50 on VoronoiDelaunay.jl which handles the space
# mapping automatically.

"""
    fcn2geom(z, rmin, rmax, imin, imax)

Linearly scale complex `z` within real range (`rmin`, `rmax`) and imaginary range (`imin`,
`imax`) into `RootsAndPoles` domain.
"""
function fcn2geom(z, rmin, rmax, imin, imax)
    zr = fcn2geom(real(z), rmin, rmax)
    zi = fcn2geom(imag(z), imin, imax)
    return complex(zr, zi)
end

"""
    geom2fcn(g, xmin, xmax)

Linearly scale `g`, which exists in the `RootsAndPoles` domain into the range of the
original function, bounded by (`xmin`, `xmax`).
"""
function geom2fcn(g, xmin, xmax)
    return (xmax - xmin)*(g - MINCOORD)/(MAXCOORD - MINCOORD) + xmin
end

"""
    geom2fcn(z, rmin, rmax, imin, imax)

Scale complex value `z` from the `RootsAndPoles` geometry domain to the comlex function""
domain with real range (`rmin`, `rmax`) and imaginary range (`imin`, `imax`).
"""
geom2fcn(z, rmin, rmax, imin, imax) = geom2fcn(real(z), imag(z), rmin, rmax, imin, imax)

"""
    geom2fcn(zr, zi, rmin, rmax, imin, imax)

Scale real and complex components `zr`, `zi` from the `RootsAndPoles` geometry domain to the
function domain.
"""
geom2fcn(zr, zi, rmin, rmax, imin, imax) = complex(geom2fcn(zr, rmin, rmax),
                                                   geom2fcn(zi, imin, imax))

"""
    geom2fcn(pt::AbstractPoint2D, rmin, rmax, imin, imax)

Scale `pt` from the geometry domain to the function domain.
"""
geom2fcn(pt::AbstractPoint2D, rmin, rmax, imin, imax) = geom2fcn(getx(pt), gety(pt), rmin,
                                                                 rmax, imin, imax)

"""
    geom2fcn(edge::DelaunayEdge, rmin, rmax, imin, imax)

Scale both points of `edge` from the geometry domain to the function domain.
"""
function geom2fcn(edge::DelaunayEdge{P}, rmin, rmax, imin, imax) where P <: AbstractPoint2D
    return (geom2fcn(geta(edge), rmin, rmax, imin, imax),
        geom2fcn(getb(edge), rmin, rmax, imin, imax))
end
