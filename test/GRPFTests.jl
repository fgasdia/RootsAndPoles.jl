using Test
using Profile
# using Gadfly
using PlotlyJS
using SpecialFunctions

include("../src/GRPF.jl")

"""
Linearly map function values within domain from `min_coord` to `max_coord`.
"""
function mapfunctionval(z, ra, rb, ia, ib)
    zr = ra*real(z) + rb
    zi = ia*imag(z) + ib
    complex(zr, zi)
end
function mapfunctionval!(z, ra, rb, ia, ib)
    for ii in eachindex(z)
        z[ii] = mapfunctionval(z[ii], ra, rb, ia, ib)
    end
end

"""
Linearly map geometry values ∈ {`min_coord`, `max_coord`} to domain bounds.

Also, there are floating point imprecisions when converting back and forth.
"""
function geom2fcn(pt::AbstractPoint2D, ra, rb, ia, ib)
    complex((getx(pt) - rb)/ra, (gety(pt) - ib)/ia)
end
geom2fcn(edge::DelaunayEdge, ra, rb, ia, ib) = (geom2fcn(geta(edge), ra, rb, ia, ib), geom2fcn(getb(edge), ra, rb, ia, ib))


@testset "Simple Rational Function" begin
    function simplefcn(z)
        w = (z - 1)*(z - im)^2*(z + 1)^3/(z + im)
    end

    # Analysis parameters
    xb = -2.  # real part begin
    xe = 2.  # real part end
    yb = -2.  # imag part begin
    ye = 2.  # imag part end
    r = 0.1  # initial mesh step
    tolerance = 1e-9

    origcoords = rectangulardomain(complex(xb, yb), complex(xe, ye), r)

    rmin, rmax = minimum(real(origcoords)), maximum(real(origcoords))
    imin, imax = minimum(imag(origcoords)), maximum(imag(origcoords))

    ra = (max_coord-min_coord)/(rmax-rmin)
    rb = max_coord - ra*rmax

    ia = (max_coord-min_coord)/(imax-imin)
    ib = max_coord - ia*imax

    mapfunctionval!(origcoords, ra, rb, ia, ib)
    newnodes = [IndexablePoint2D(real(coord), imag(coord), idx) for (idx, coord) in enumerate(origcoords)]
    tess = DelaunayTessellation2D{IndexablePoint2D}(2000)

    tess, 𝓔, quadrants = tesselate!(tess, newnodes, pt -> simplefcn(geom2fcn(pt, ra, rb, ia, ib)),
                                    e -> geom2fcn(e, ra, rb, ia, ib), tolerance)

    # x, y = getplotxy(delaunayedges(tess))
    # t1 = PlotlyJS.scatter(;x=x, y=y) # mode="markers", marker_size=3)
    # l = Layout(width=600, height=600)
    # PlotlyJS.plot(t1, l)

    # rz = geom2fcn.(regions[1])
    # t1 = PlotlyJS.scatter(;x=real(rz), y=imag(rz), mode="markers", marker_size=6)
    # l = Layout(width=600, height=600)
    # PlotlyJS.plot(t1, l)

    𝐶 = contouredges(tess, 𝓔)
    regions = evaluateregions!(𝐶, e -> geom2fcn(e, ra, rb, ia, ib))

    zroots, zroots_multiplicity, zpoles, zpoles_multiplicity = rootsandpoles(regions, quadrants, e -> geom2fcn(e, ra, rb, ia, ib))
    println("zroots: ")
    display(zroots)
    println("\nzpoles: ")
    display(zpoles)
    println()

    @test count(real(zroots) .> 0) == 2
    @test count(imag(zroots) .> 0) == 2
end

@testset "Complex Modes" begin
    function complexmodes(z)
        z *= 10
        f = 5e9
        c = 3e8
        μ₀ = 4e-7π
        ϵ₀ = 1e-9/36/π
        a = 6.35e-3
        b = 10e-3
        ϵᵣ₁ = 10
        ϵᵣ₂ = 1
        m = 1

        ω = 2π*f
        k₀ = ω/c
        α = real(z)*k₀
        β = imag(z)*k₀
        γ = α + im*β
        ϵ₁ = ϵ₀*ϵᵣ₁
        ϵ₂ = ϵ₀*ϵᵣ₂
        μ₁ = μ₀
        μ₂ = μ₀
        κ₁ = sqrt(γ^2 + k₀^2*ϵᵣ₁)
        κ₂ = sqrt(γ^2 + k₀^2*ϵᵣ₂)
        η₁ = sqrt(μ₁/ϵ₁)
        η₂ = sqrt(μ₂/ϵ₂)

        Jm_a1 = besselj1(κ₁*a)
        Jm_a2 = besselj1(κ₂*a)
        Ym_a2 = bessely1(κ₂*a)
        Jm_b2 = besselj1(κ₂*b)
        Ym_b2 = bessely1(κ₂*b)
        DJm_a1 = (besselj0(κ₁*a) - besselj(2, κ₁*a))/2
        DJm_a2 = (besselj0(κ₂*a) - besselj(2, κ₂*a))/2
        DJm_b2 = (besselj0(κ₂*b) - besselj(2, κ₂*b))/2
        DYm_a2 = (bessely0(κ₂*a) - bessely(2, κ₂*a))/2
        DYm_b2 = (bessely0(κ₂*b) - bessely(2, κ₂*b))/2

        W = @SMatrix [Jm_a1                 0                       -Jm_a2              -Ym_a2              0                       0;
                      0                     Jm_a1/η₁                0                   0                   -Jm_a2/η₂               -Ym_a2/η₂;
                      γ*m*Jm_a1/(a*κ₁^2)    -ω*μ₁*DJm_a1/(κ₁*η₁)    -γ*m*Jm_a2/(a*κ₂^2) -γ*m*Ym_a2/(a*κ₂^2) ω*μ₂*DJm_a2/(κ₂*η₂)     ω*μ₂*DYm_a2/(κ₂*η₂);
                      -ω*ϵ₁*DJm_a1/κ₁       -m*γ*Jm_a1/(a*κ₁^2*η₁)  ω*ϵ₂*DJm_a2/κ₂      ω*ϵ₂*DYm_a2/κ₂      m*γ*Jm_a2/(a*κ₂^2*η₂)   m*γ*Ym_a2/(a*κ₂^2*η₂);
                      0                     0                       Jm_b2               Ym_b2               0                       0;
                      0                     0                       γ*m*Jm_b2/(b*κ₂^2)  γ*m*Ym_b2/(b*κ₂^2)  -ω*μ₂*DJm_b2/(κ₂*η₂)    -ω*μ₂*DYm_b2/(κ₂*η₂)]
        w = det(W)
    end

    R = 1.
    r = 0.15

    origcoords = diskdomain(R, r)
    tolerance = 1e-9

    rmin, rmax = minimum(real(origcoords)), maximum(real(origcoords))
    imin, imax = minimum(imag(origcoords)), maximum(imag(origcoords))

    ra = (max_coord-min_coord)/(rmax-rmin)
    rb = max_coord - ra*rmax

    ia = (max_coord-min_coord)/(imax-imin)
    ib = max_coord - ia*imax

    mapfunctionval!(origcoords, ra, rb, ia, ib)
    newnodes = [IndexablePoint2D(real(coord), imag(coord), idx) for (idx, coord) in enumerate(origcoords)]
    tess = DelaunayTessellation2D{IndexablePoint2D}(2000)

    tess, 𝓔, quadrants = tesselate!(tess, newnodes, pt -> complexmodes(geom2fcn(pt, ra, rb, ia, ib)),
                                    e -> geom2fcn(e, ra, rb, ia, ib), tolerance)

    𝐶 = contouredges(tess, 𝓔)
    regions = evaluateregions!(𝐶, e -> geom2fcn(e, ra, rb, ia, ib))

    zroots, zroots_multiplicity, zpoles, zpoles_multiplicity = rootsandpoles(regions, quadrants, e -> geom2fcn(e, ra, rb, ia, ib))
    println("zroots: ")
    display(zroots)
    println("\nzpoles: ")
    display(zpoles)
    println()

    @test count(real(zroots) .> 0) == 6
    @test count(imag(zroots) .> 0) == 5
end

@testset "Lossy Multilayered Waveguide" begin
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
    xb = 1.  # real part begin
    xe = 2.5  # real part end
    yb = -1.  # imag part begin
    ye = 1.  # imag part end
    r = 0.5  # initial mesh step
    tolerance = 1e-9

    origcoords = rectangulardomain(complex(xb, yb), complex(xe, ye), r)

    rmin, rmax = minimum(real(origcoords)), maximum(real(origcoords))
    imin, imax = minimum(imag(origcoords)), maximum(imag(origcoords))

    ra = (max_coord-min_coord)/(rmax-rmin)
    rb = max_coord - ra*rmax

    ia = (max_coord-min_coord)/(imax-imin)
    ib = max_coord - ia*imax

    mapfunctionval!(origcoords, ra, rb, ia, ib)
    newnodes = [IndexablePoint2D(real(coord), imag(coord), idx) for (idx, coord) in enumerate(origcoords)]
    tess = DelaunayTessellation2D{IndexablePoint2D}(2000)

    tess, 𝓔, quadrants = tesselate!(tess, newnodes, pt -> wvgd(geom2fcn(pt, ra, rb, ia, ib)),
                                    e -> geom2fcn(e, ra, rb, ia, ib), tolerance)

    # x, y = getplotxy(delaunayedges(tess))
    # t1 = PlotlyJS.scatter(;x=x, y=y) # mode="markers", marker_size=3)
    # l = Layout(width=600, height=600)
    # PlotlyJS.plot(t1, l)
    #
    # rz = geom2fcn.(regions[1])
    # t1 = PlotlyJS.scatter(;x=real(rz), y=imag(rz), mode="markers", marker_size=6)
    # l = Layout(width=600, height=600)
    # PlotlyJS.plot(t1, l)

    𝐶 = contouredges(tess, 𝓔)
    regions = evaluateregions!(𝐶, e -> geom2fcn(e, ra, rb, ia, ib))

    zroots, zroots_multiplicity, zpoles, zpoles_multiplicity = rootsandpoles(regions, quadrants, e -> geom2fcn(e, ra, rb, ia, ib))
    println("zroots: ")
    display(zroots)
    println("\nzpoles: ")
    display(zpoles)
    println()

    @test count(real(zroots) .> 0) == 7
    @test count(imag(zroots) .> 0) == 0
end

@testset "Graphene Transmission Line" begin
    function graphenefunction(z)
        f = 1e12
        c = 299792458.
        μ₀ = 4π*1e-7
        ϵ₀ = 1/(μ₀*c^2)

        e = 1.602176565e-19
        kB = 1.3806488e-23
        hk = 1.05457168e-34
        vFe = 1e6
        muc = 0.05*e
        t = 0.135e-12
        T = 300
        ϵᵣ₁ = 1.
        ϵᵣ₂ = 11.9

        ω = 2π*f
        k₀ = ω/c
        kᵣ₀ = -im*z*k₀

        Slo=-im*e^2*kB*T*log(2+2*cosh(muc/kB/T)) / (π*hk^2*(ω-im/t))

        a = -3*vFe^2*Slo/(4*(ω-im/t)^2)
        b = a/3

        Y1TM = ω*ϵᵣ₁*ϵ₀/sqrt(ϵᵣ₁*k₀^2 - kᵣ₀^2);
        Y2TM = ω*ϵᵣ₂*ϵ₀/sqrt(ϵᵣ₂*k₀^2 - kᵣ₀^2);
        YSTM = Slo + 1*a*kᵣ₀^2 + 1*b*kᵣ₀^2;

        w = (Y1TM + Y2TM + YSTM)*(-Y1TM + Y2TM + YSTM)*(Y1TM - Y2TM + YSTM)*(-Y1TM - Y2TM + YSTM) # four Riemann sheets
    end

    # Analysis parameters
    xb = -100.  # real part begin
    xe = 400.  # real part end
    yb = -100.  # imag part begin
    ye = 400.  # imag part end
    r = 18.  # initial mesh step
    tolerance = 1e-9

    origcoords = rectangulardomain(complex(xb, yb), complex(xe, ye), r)

    rmin, rmax = minimum(real(origcoords)), maximum(real(origcoords))
    imin, imax = minimum(imag(origcoords)), maximum(imag(origcoords))

    ra = (max_coord-min_coord)/(rmax-rmin)
    rb = max_coord - ra*rmax

    ia = (max_coord-min_coord)/(imax-imin)
    ib = max_coord - ia*imax

    origcoords = mapfunctionval.(origcoords, ra, rb, ia, ib)
    newnodes = [IndexablePoint2D(real(coord), imag(coord), idx) for (idx, coord) in enumerate(origcoords)]
    tess = DelaunayTessellation2D{IndexablePoint2D}(2000)

    tess, 𝓔, quadrants = tesselate!(tess, newnodes, pt -> graphenefunction(geom2fcn(pt, ra, rb, ia, ib)),
                                    e -> geom2fcn(e, ra, rb, ia, ib), tolerance)

#
#     # x, y = getplotxy(delaunayedges(tess))
#     # set_default_plot_size(15cm, 15cm)
#     # p = plot(x=x, y=y, Geom.path, Scale.x_continuous(minvalue=1.0, maxvalue=2.0), Scale.y_continuous(minvalue=1.0, maxvalue=2.0))
#     # draw(SVG("graphenefunction.svg", 6inch, 6inch), p)
#
#     𝐶 = contouredges(tess, 𝓔)
#     regions = evaluateregions!(𝐶, e -> geom2fcn(e, ra, rb, ia, ib))
#
#     zroots, zroots_multiplicity, zpoles, zpoles_multiplicity = rootsandpoles(regions, quadrants, e -> geom2fcn(e, ra, rb, ia, ib))
#     println("zroots: ")
#     display(zroots)
#     println("\nzpoles: ")
#     display(zpoles)
#     println()
#
#     @test count(real(zroots) .> 0) == 6
#     @test count(imag(zroots) .> 0) == 6
end

@testset "Default" begin
    function defaultfcn(z)
        f = 1e9
        ϵᵣ = 5 - 2im
        μᵣ = 1 - 2im
        d = 1e-2
        c = 3e8
        ω = 2π*f
        k₀ = ω/c
        cc = ϵᵣ^2*(k₀*d)^2*(ϵᵣ*μᵣ - 1)
        w = ϵᵣ^2*z^2 + z^2*tan(z)^2 - cc
    end

    # Analysis parameters
    xb = -2.  # real part begin
    xe = 2.  # real part end
    yb = -2.  # imag part begin
    ye = 2.  # imag part end
    r = 0.2  # initial mesh step
    tolerance = 1e-9

    origcoords = rectangulardomain(complex(xb, yb), complex(xe, ye), r)

    rmin, rmax = minimum(real(origcoords)), maximum(real(origcoords))
    imin, imax = minimum(imag(origcoords)), maximum(imag(origcoords))

    ra = (max_coord-min_coord)/(rmax-rmin)
    rb = max_coord - ra*rmax

    ia = (max_coord-min_coord)/(imax-imin)
    ib = max_coord - ia*imax

    mapfunctionval!(origcoords, ra, rb, ia, ib)
    newnodes = [IndexablePoint2D(real(coord), imag(coord), idx) for (idx, coord) in enumerate(origcoords)]
    tess = DelaunayTessellation2D{IndexablePoint2D}(2000)

    @time begin
        tess, 𝓔, quadrants = tesselate!(tess, newnodes, pt -> defaultfcn(geom2fcn(pt, ra, rb, ia, ib)),
                                        e -> geom2fcn(e, ra, rb, ia, ib), tolerance)

        𝐶 = contouredges(tess, 𝓔)
        regions = evaluateregions!(𝐶, e -> geom2fcn(e, ra, rb, ia, ib))

        zroots, zroots_multiplicity, zpoles, zpoles_multiplicity = rootsandpoles(regions, quadrants, e -> geom2fcn(e, ra, rb, ia, ib))
    end
    println("zroots: ")
    display(zroots)
    println("\nzpoles: ")
    display(zpoles)
    println()

    @test count(real(zroots) .> 0) == 3
    @test count(imag(zroots) .> 0) == 3
end
