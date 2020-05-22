function distance(p1::Complex, p2::Complex)
    return sqrt((real(p2) - real(p1))^2 + (imag(p2) - imag(p1))^2)
end
distance((p1, p2)::Tuple{Complex, Complex}) = distance(p1, p2)
function distance(p1::AbstractPoint2D, p2::AbstractPoint2D)
    return sqrt((getx(p2)-getx(p1))^2 + (gety(p2)-gety(p1))^2)
end
distance(e::DelaunayEdge) = distance(getb(e), geta(e))
function distance(p1arr::AbstractArray{AbstractPoint2D}, p2::AbstractPoint2D)
    return sqrt.((getx.(p1arr).-getx(p2)).^2 .+ (gety.(p1arr).-gety(p2)).^2)
end

"""
    longedge(edge, tolerance, geom2fcn)

Return true if `edge` has length greater than `tolerance`.
"""
function longedge(edge::DelaunayEdge, tolerance, g2f::Geometry2Function)
    return distance(g2f(edge)) > tolerance
end

"""
    fcn2geom(z, ra, rb, ia, ib)

Linearly map function values `z` within domain from `ra` to `rb` and `ia` to `ib`. Necessary
because [VoronoiDelaunay.jl](https://github.com/JuliaGeometry/VoronoiDelaunay.jl) requires
point coordinates be within `min_coord <= x <= max_coord` where `min_coord=1.0+eps(Float64)`
and `max_coord=2.0-2eps(Float64)`. `max_coord` and `min_coord` are provided by
`VoronoiDelaunay.jl`
"""
function fcn2geom(z, ra, rb, ia, ib)
    zr = ra*real(z) + rb
    zi = ia*imag(z) + ib
    return complex(zr, zi)
end

"""
    geom2fcn(pt, ra, rb, ia, ib)

Linearly map geometry values ∈ {`min_coord`, `max_coord`} to domain bounds.

Note: There are floating point errors when converting back and forth.
"""
function geom2fcn(pt::AbstractPoint2D, ra, rb, ia, ib)
    return complex((getx(pt) - rb)/ra, (gety(pt) - ib)/ia)
end
function geom2fcn(edge::VoronoiDelaunay.DelaunayEdge{P}, ra, rb, ia, ib) where P <: AbstractPoint2D
    return geom2fcn(geta(edge), ra, rb, ia, ib), geom2fcn(getb(edge), ra, rb, ia, ib)
end
function geom2fcn(x, y, ra, rb, ia, ib)
    return complex((x-rb)/ra, (y-ib)/ia)
end
