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
    rectangulardomain(Zb, Ze, Δr)

Generate initial mesh node coordinates for a rectangular domain ∈ {`Zb` to `Ze`} with
initial mesh step `Δr`.

See also: `rect_dom.m`
"""
function rectangulardomain(Zb::Complex, Ze::Complex, Δr)
    X = real(Ze) - real(Zb)
    Y = imag(Ze) - imag(Zb)

    # BUG?: Why are n and m in real and imag (x and y) different?
    n = ceil(Int, Y/Δr + 1)
    dy = Y/(n-1)
    m = ceil(Int, X/sqrt(Δr^2 - dy^2/4) + 1)
    dx = X/(m-1)

    vx = range(real(Zb), length=m, stop=real(Ze))
    vy = range(imag(Zb), length=n, stop=imag(Ze))

    tmp = ones(n)
    tmp[n] = 0.0

    # TODO: This looks _so_ matlaby
    y = repeat(vy, m)
    y .+= 0.5*dy*kron((1 .+ (-1).^(1:m))/2, tmp)
    y = [y; fill(imag(Zb), size(2:2:m))]
    x = reshape(repeat(vx', n), m*n)
    x = [x; ((2:2:m) .- 1)*dx .+ real(Zb)]

    # NOTE: Matlab values of `x` differ slightly because of Matlab's float handling and
    # transpose operator.
    # `sum(x)` is much closer between Julia and Matlab if the above line for `x` is:
    # x = [x' ((2:2:m) .- 1)'*dx .+ real(Zb)]'
    # x = reshape(x, length(x))

    return complex.(x, y)
end

"""
    diskdomain(R, Δr)

Generate initial mesh coordinates for a circular disk domain of radius `R` and center (0, 0)
for ``|z| < R``.

See also: `disk_dom.m`
"""
function diskdomain(R, Δr)
    h = Δr*sqrt(3)/2
    n = 1 + round(Int, R/h)
    Rn = (1:n)*R/n
    newnodes = [complex(0.0)]
    f₀ = 0.0
    np = 6
    for ii = 1:n
        f = f₀ .+ range(0, stop=2π, length=np+1)
        xyn = Rn[ii]*cis.(f[1:end-1])
        append!(newnodes, xyn)
        f₀ += π/6/n
        np += 6
    end
    return newnodes
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
