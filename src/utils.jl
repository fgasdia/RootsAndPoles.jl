distance(p1::Complex, p2::Complex) = hypot(p2 - p1)
distance((p1, p2)::Tuple{Complex, Complex}) = distance(p1, p2)
function distance(p1::AbstractPoint2D, p2::AbstractPoint2D)
    z1 = complex(p1.x, p1.y)
    z2 = complex(p2.x, p2.y)
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

"""
    getplotdata(tess, quadrants, phasediffs, g2f)

Return complex coordinates of every tesselation edge and a vector with the corresponding
index specifiying its phase quadrant.

This function is useful for plotting the final tesselation using results from calling
`grpf` with the `PlotData()` argument.

The edge colors are (where Q means "quadrant"):

| Color index |         Meaning         |  Quadrant phase   |
|:-----------:|:-----------------------:|:-----------------:|
|      1      |           Q 1           | 0 ≤ arg f < π/2   |
|      2      |           Q 2           | π/2 ≤ arg f < π   |
|      3      |           Q 3           | π ≤ arg f < 3π/2  |
|      4      |           Q 4           | 3π/2 ≤ arg f < 2π |
|      5      | phase change (boundary) |         -         |
|      6      | candidate edge          |         -         |

# Example

Although not self contained, this example shows how `getplotdata` can be used:

```julia
roots, poles, quads, phasediffs, tess, g2f = grpf(f, origcoords, PlotData());
z, edgecolors = getplotdata(tess, quadrants, phasediffs, g2f);

using Plots
plot(real(z), imag(z), group=edgecolors)
```
"""
function getplotdata(tess, quadrants, phasediffs, g2f)
    edges = [e for e in delaunayedges(tess)]

    edgecolors = fill(5, length(edges))
    for ii in eachindex(edges)
        if phasediffs[ii] == 2  # candidate edges
            edgecolors[ii] = 6
        elseif phasediffs[ii] == 0
            edgecolors[ii] = quadrants[getindex(geta(edges[ii]))]
        end
    end

    edgecolors = repeat(edgecolors, inner=3)  # x, y is xₐ, xᵦ, NaN, repeat

    x, y = getplotxy(edges)
    z = g2f.(x, y)

    I = sortperm(edgecolors)

    return z[I], edgecolors[I]
end
