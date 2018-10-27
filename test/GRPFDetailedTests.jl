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


@testset "Default function" begin
    function testfunction(z)
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
    testfunction(z::AbstractArray) = [testfunction(zz) for zz in z]


    # Analysis parameters
    xb = -2.  # real part begin
    xe = 2.  # real part end
    yb = -2.  # imag part begin
    ye = 2.  # imag part end
    r = 0.2  # initial mesh step

    origcoords = rectangulardomain(complex(xb, yb), complex(xe, ye), r)

    tolerance = 1e-9
    skinnytriangle = 3

    rmin, rmax = minimum(real(origcoords)), maximum(real(origcoords))
    imin, imax = minimum(imag(origcoords)), maximum(imag(origcoords))

    ra = (max_coord-min_coord)/(rmax-rmin)
    rb = max_coord - ra*rmax

    ia = (max_coord-min_coord)/(imax-imin)
    ib = max_coord - ia*imax

    mapfunctionval!(origcoords, ra, rb, ia, ib)
    newnodes = [IndexablePoint2D(real(coord), imag(coord), idx) for (idx, coord) in enumerate(origcoords)]
    tess = DelaunayTessellation2D{IndexablePoint2D}(2000)

    quadrants = Vector{Int64}(undef, length(newnodes))
    # quadrants = Array{Int64}(undef, 0, 2)
    assignquadrants!(quadrants, newnodes, pt -> testfunction(geom2fcn(pt, ra, rb, ia, ib)))

    @test count(quadrants .== 1) == 221
    @test count(quadrants .== 2) == 106
    @test count(quadrants .== 3) == 124
    @test count(quadrants .== 4) == 86

    push!(tess, newnodes)
    𝓔, phasediffs = candidateedges(tess, quadrants)
    @test length(phasediffs) == 1520
    @test count(phasediffs .== 0) == 1330
    @test count(phasediffs .== 1) == 92
    @test count(phasediffs .== 2) == 10
    @test count(phasediffs .== 3) == 88

    matlab𝓔 = [52 73;
               53 73;
               72 73;
               73 95;
               182 202;
               323 344;
               431 452;
               452 453;
               452 473;
               452 474]

    # Sorting for Matlab comparison only
    testedges = sort([getindex.(geta.(𝓔)) getindex.(getb.(𝓔))], dims=2)
    testedges = sort(testedges, dims=1, by=x->x[1])

    @test matlab𝓔 == testedges

    select𝓔 = selectedges(𝓔, tolerance, e -> geom2fcn(e, ra, rb, ia, ib))
    newE = filter(e -> distance(geom2fcn(e, ra, rb, ia, ib)...) > tolerance, 𝓔)

    tttfcn(e) = distance(geom2fcn(e, ra, rb, ia, ib)...) > tolerance
    select𝓔 = filter(tttfcn, 𝓔)

    testedges = sort([getindex.(geta.(select𝓔)) getindex.(getb.(select𝓔))], dims=2)
    testedges = sort(testedges, dims=1, by=x->x[1])

    @test matlab𝓔 == testedges  # In first loop, they're the same
    @test minimum(select𝓔) ≈ 0.19436506316151
    @test maximum(select𝓔) ≈ 0.20000000000000

    trianglecounts = counttriangleswithnodes(tess, select𝓔)

    @test count(trianglecounts .> 1) == 20
    @test count(trianglecounts .== 1) == 40
    @test maximum(trianglecounts) == 3
    @test sum(trianglecounts) == 84

    zone1triangles, zone2triangles = splittriangles(tess, trianglecounts)

    # zone1triangles = [tr for (idx, tr) in enumerate(tess) if trianglecounts[idx] > 1]
    # zone1count = length(zone1triangles)
    # @test zone1count == 20

    newnodes = Vector{IndexablePoint2D}()
    zone1newnodes!(newnodes, zone1triangles, e -> geom2fcn(e, ra, rb, ia, ib), tolerance)
    @test length(newnodes) == 42

    # zone2triangles = [tr for (idx, tr) in enumerate(tess) if trianglecounts[idx] == 1]
    # @test length(zone2triangles) == 40

    zone2newnodes!(newnodes, zone2triangles)
    @test length(newnodes) == 42

    𝐶 = contouredges(tess, 𝓔)

    regions = evaluateregions!(𝐶, geom2fcn)

    zroots, zroots_multiplicity, zpoles, zpoles_multiplicity = rootsandpoles(regions, quadrants, e -> geom2fcn(e, ra, rb, ia, ib))
end


@testset "Simple rational function breakdown" begin
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

    # Initialize
    numnodes = tess._total_points_added
    @assert numnodes == 0

    𝓔 = Vector{DelaunayEdge}()
    quadrants = Vector{Int64}()

    iteration = 0
    while (iteration < maxiterations) & (numnodes < maxnodes)
        iteration += 1
        println("Iteration: ", iteration)

        # Determine which quadrant function value belongs at each node
        numnewnodes = length(newnodes)
        append!(quadrants, Vector{Int64}(undef, numnewnodes))
        # if iteration > 1
        #     nn = [t for t in tess if true]
        #     nna = geta.(nn)
        #     nnb = getb.(nn)
        #     nnc = getc.(nn)
        #     nnn = unique([nna; nnb; nnc])
        #     @assert length(quadrants) == maximum(getindex.(nnn))
        #     @assert allunique(getindex.(nnn))
        # end
        assignquadrants!(quadrants, newnodes, pt -> simplefcn(geom2fcn(pt, ra, rb, ia, ib)))
        if iteration == 1
            @test sum(quadrants) == 5020
            @test count(quadrants .== 1) == 509
            @test count(quadrants .== 2) == 495
            @test count(quadrants .== 3) == 431
            @test count(quadrants .== 4) == 557
        elseif iteration == 2
            @test sum(quadrants) == 5148
            @test count(quadrants .== 1) == 522
            @test count(quadrants .== 2) == 507
            @test count(quadrants .== 3) == 444
            @test count(quadrants .== 4) == 570
        elseif iteration == 3
            @test sum(quadrants) == 5314
            @test count(quadrants .== 1) == 537
            @test count(quadrants .== 2) == 527
            @test count(quadrants .== 3) == 465
            @test count(quadrants .== 4) == 582
        end

        # t1 = PlotlyJS.scatter(;x=getx.(newnodes), y=gety.(newnodes), mode="markers", marker_size=3)
        # l = Layout(width=600, height=600)
        # PlotlyJS.plot(t1, l)

        # Add new nodes to `tess`
        push!(tess, newnodes)
        numnodes += numnewnodes
        iteration == 1 && @test numnodes == 1992
        iteration == 2 && @test numnodes == 2043
        iteration == 3 && @test numnodes == 2111

        # Determine candidate edges that may be near a root or pole
        𝓔, phasediffs = candidateedges(tess, quadrants)
        isempty(𝓔) && error("No roots in the domain")
        if iteration == 1
            # @test sum(phasediffs) == 1927
            # @test count(phasediffs .== 0) == 4916
            # @test count(phasediffs .== 1) == 354
            # @test count(phasediffs .== 2) == 11
            # @test count(phasediffs .== 3) == 517
            @test length(𝓔) == 11
        elseif iteration == 2
            # @test sum(phasediffs) == 2054
            # @test count(phasediffs .== 0) == 4998
            # @test count(phasediffs .== 1) == 395
            # @test count(phasediffs .== 2) == 15
            # @test count(phasediffs .== 3) == 543
            @test length(𝓔) == 15
        elseif iteration == 3
            @test length(𝓔) == 12
            # @test sum(phasediffs) == 2247
        end

        # Select candidate edges that are longer than the chosen tolerance
        select𝓔, min𝓔length, max𝓔length = selectedges(𝓔, tolerance, e -> geom2fcn(e, ra, rb, ia, ib))
        @debug "Candidate edges length min: $min𝓔length, max: $max𝓔length"
        max𝓔length < tolerance && return tess, 𝓔, quadrants
        if iteration == 1
            @test min𝓔length ≈ 0.0987071244830946 atol=1e-14
            @test max𝓔length ≈ 0.1000000000000001 atol=1e-14
            @test select𝓔 == 𝓔
        elseif iteration == 2
            @test min𝓔length ≈ 0.0493535622415472 atol=1e-14
            @test max𝓔length ≈ 0.0500000000000000 atol=1e-14
            @test select𝓔 == 𝓔
        elseif iteration == 3
            @test min𝓔length ≈ 0.0246767811207736 atol=1e-14
            @test max𝓔length ≈ 0.0250000000000000 atol=1e-14
        end

        # How many times does each triangle contain a `select𝓔` node?
        trianglecounts = counttriangleswithnodes(tess, select𝓔)
        if iteration == 1
            @test length(trianglecounts) == 3807
            @test sum(trianglecounts) == 90
            @test count(trianglecounts .== 1) == 39
            @test count(trianglecounts .== 2) == 21
            @test count(trianglecounts .== 3) == 3
        elseif iteration == 2
            @test length(trianglecounts) == 3909
            @test sum(trianglecounts) == 123
            @test count(trianglecounts .== 1) == 40
            @test count(trianglecounts .== 2) == 19
            @test count(trianglecounts .== 3) == 15
        end
        zone1triangles, zone2triangles = splittriangles(tess, trianglecounts)
        if iteration == 1
            @test length(zone1triangles) == 24
            @test length(zone2triangles) == 39
        elseif iteration == 2
            @test length(zone1triangles) == 34
            @test length(zone2triangles) == 40
        end

        # Add new nodes in zone 1
        newnodes = Vector{IndexablePoint2D}()
        zone1newnodes!(newnodes, zone1triangles, e -> geom2fcn(e, ra, rb, ia, ib), tolerance)
        testnodes1 = sort(geom2fcn.(newnodes, ra, rb, ia, ib), by=x->real(x))
        if iteration == 1
            @test length(newnodes) == 51
            @test sum(testnodes1) ≈ complex(-19.021276595744673, 12.450000000000003) atol=1e-13
            @test count(real.(testnodes1) .> 0) == 10
            @test count(imag.(testnodes1) .> 0) == 30
            @test count(real.(testnodes1) .<= 0) == 41  # XXX: Matlab doesn't have any 0s, they are very small negative numbers. does this change results elsewhere?
            @test count(imag.(testnodes1) .<= 0) == 21
        elseif iteration == 2
            @test length(newnodes) == 68
            @test sum(testnodes1) ≈ complex(-33.021276595744681, 15.550000000000001) atol=1e-13
            @test count(real.(testnodes1) .> 0) == 14
            @test count(imag.(testnodes1) .> 0) == 42
            @test count(real.(testnodes1) .<= 0) == 54  # XXX: Matlab doesn't have any 0s, they are very small negative numbers. does this change results elsewhere?
            @test count(imag.(testnodes1) .<= 0) == 26
        end

        # Add new nodes in zone 2
        zone2newnodes!(newnodes, zone2triangles)
        testnodes2 = sort(geom2fcn.(newnodes, ra, rb, ia, ib), by=x->real(x))
        if iteration == 1
            @test length(newnodes) == 51
            @test sum(testnodes2) ≈ complex(-19.021276595744673, 12.450000000000003)
            @test count(real.(testnodes2) .> 0) == 10
            @test count(imag.(testnodes2) .> 0) == 30
            @test count(real.(testnodes2) .<= 0) == 41
            @test count(imag.(testnodes2) .<= 0) == 21
        elseif iteration == 2
            @test length(newnodes) == 68
            @test sum(testnodes2) ≈ complex(-33.02127659574468, 15.550000000000001)
            @test count(real.(testnodes2) .> 0) == 14
            @test count(imag.(testnodes2) .> 0) == 42
            @test count(real.(testnodes2) .<= 0) == 54
            @test count(imag.(testnodes2) .<= 0) == 26
        end

        # Have to assign indexes to new nodes (which are all currently -1)
        setindex!.(newnodes, (1:length(newnodes)).+numnodes)
    end
end


open("newnodes2.csv", "w") do f
    for ii in eachindex(newnodes)
        z = geom2fcn(newnodes[ii], ra, rb, ia, ib)
        write(f, string(real(z)), ", ", string(imag(z)), "\n")
    end
end

open("fvals2.csv", "w") do f
    for ii in eachindex(newnodes)
        z = simplefcn(geom2fcn(newnodes[ii], ra, rb, ia, ib))
        write(f, string(real(z)), ", ", string(imag(z)), "\n")
    end
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

    R = 1
    r = 0.15

    origcoords = diskdomain(R, r)
    tolerance = 1e-9

    @test length(origcoords) == 271
    @test maximum(real(origcoords)) == 1
    @test maximum(imag(origcoords)) ≈ 0.998308158271268 atol=1e-15
    @test minimum(real(origcoords)) == -1
    @test minimum(imag(origcoords)) ≈ -0.998308158271268 atol=1e-15
    @test sum(real(origcoords)) ≈ -1.243449787580175e-14 atol=1e-14
    @test sum(imag(origcoords)) ≈ -4.829470157119431e-15 atol=1e-14

    rmin, rmax = minimum(real(origcoords)), maximum(real(origcoords))
    imin, imax = minimum(imag(origcoords)), maximum(imag(origcoords))

    ra = (max_coord-min_coord)/(rmax-rmin)
    rb = max_coord - ra*rmax

    ia = (max_coord-min_coord)/(imax-imin)
    ib = max_coord - ia*imax

    mapfunctionval!(origcoords, ra, rb, ia, ib)
    newnodes = [IndexablePoint2D(real(coord), imag(coord), idx) for (idx, coord) in enumerate(origcoords)]
    tess = DelaunayTessellation2D{IndexablePoint2D}(2000)

    #####

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
        assignquadrants!(quadrants, newnodes, pt -> complexmodes(geom2fcn(pt, ra, rb, ia, ib)))

        if iteration == 1
            @test length(quadrants) == 271
            @test count(quadrants .== 1) == 60
            @test count(quadrants .== 2) == 76
            @test count(quadrants .== 3) == 71
            @test count(quadrants .== 4) == 64
        end

        # Add new nodes to `tess`
        push!(tess, newnodes)
        numnodes += numnewnodes

        # Determine candidate edges that may be near a root or pole
        𝓔, phasediffs = candidateedges(tess, quadrants)
        isempty(𝓔) && error("No roots in the domain")

        # Select candidate edges that are longer than the chosen tolerance
        select𝓔, min𝓔length, max𝓔length = selectedges(𝓔, tolerance, e -> geom2fcn(e, ra, rb, ia, ib))
        @debug "Candidate edges length min: $min𝓔length, max: $max𝓔length"
        max𝓔length < tolerance && return tess, 𝓔, quadrants

        # How many times does each triangle contain a `select𝓔` node?
        trianglecounts = counttriangleswithnodes(tess, select𝓔)
        zone1triangles, zone2triangles = splittriangles(tess, trianglecounts)

        # Add new nodes in zone 1
        newnodes = Vector{IndexablePoint2D}()
        zone1newnodes!(newnodes, zone1triangles, e -> geom2fcn(e, ra, rb, ia, ib), tolerance)

        # Add new nodes in zone 2
        zone2newnodes!(newnodes, zone2triangles)

        # Have to assign indexes to new nodes (which are all currently -1)
        setindex!.(newnodes, (1:length(newnodes)).+numnodes)
    end

    ####

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

    fcn(pt) = graphenefunction(geom2fcn(pt, ra, rb, ia, ib))
    geom2fcn(e) = geom2fcn(e, ra, rb, ia, ib)

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

    # tess, 𝓔, quadrants = tesselate!(tess, newnodes, pt -> graphenefunction(geom2fcn(pt, ra, rb, ia, ib)),
    #                                 e -> geom2fcn(e, ra, rb, ia, ib), tolerance)


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



#
#     # x, y = getplotxy(delaunayedges(tess))
#     # set_default_plot_size(15cm, 15cm)
#     # p = plot(x=x, y=y, Geom.path, Scale.x_continuous(minvalue=1.0, maxvalue=2.0), Scale.y_continuous(minvalue=1.0, maxvalue=2.0))
#     # draw(SVG("graphenefunction.svg", 6inch, 6inch), p)
#
    𝐶 = contouredges(tess, 𝓔)
    regions = evaluateregions!(𝐶, e -> geom2fcn(e, ra, rb, ia, ib))

    zroots, zroots_multiplicity, zpoles, zpoles_multiplicity = rootsandpoles(regions, quadrants, e -> geom2fcn(e, ra, rb, ia, ib))
    println("zroots: ")
    display(zroots)
    println("\nzpoles: ")
    display(zpoles)
    println()

    @test count(real(zroots) .> 0) == 6
    @test count(imag(zroots) .> 0) == 6
end
