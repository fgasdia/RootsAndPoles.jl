function test_GRPFParams()
    ga = GRPFParams()
    gb = GRPFParams(100, 500000, 3, 5000, 1e-9, false)
    gc = GRPFParams(200, 10000, 3, 10000, 1e-8, true)

    # equality
    @test ga == gb
    @test ga != gc

    # constructors
    @test GRPFParams(5000, 1e-3) == GRPFParams(100, 500000, 3, 5000, 1e-3, false)
    @test GRPFParams(5000, 1e-3, true) == GRPFParams(100, 500000, 3, 5000, 1e-3, true)

    @test_logs (:warn,
        "GRPFParams `tess_sizehint` is greater than `maxnodes`") GRPFParams(100, 10, 3,
                                                                    5000, 1e-9, false)
end

function test_functionstructs()
    rmin, rmax = -1.0, 1.0
    imin, imax = -1.0, 1.0

    g2f = RP.Geometry2Function(rmin, rmax, imin, imax)
    @test eltype(g2f) == Float64

    fval = complex(-1.0, 1.0)
    f2g = RP.fcn2geom(fval, rmin, rmax, imin, imax)

    @test reim(f2g) == (RP.MINCOORD, RP.MAXCOORD)
    @test g2f(f2g) == fval

    fval = complex(0.2, -0.1)
    f2g = RP.fcn2geom(fval, rmin, rmax, imin, imax)
    @test g2f(f2g) ≈ fval   atol=1e-15  # floating point limitation

    p = RP.Point2D(reim(f2g)...)
    @test g2f(p) ≈ fval     atol=1e-15

    fval2 = complex(0.4, 0.8)
    f2g2 = RP.fcn2geom(fval2, rmin, rmax, imin, imax)
    p2 = RP.Point2D(reim(f2g2)...)
    e = RP.DelaunayEdge(p, p2)
    f1, f2 = g2f(e)
    @test f1 ≈ fval     atol=1e-15
    @test f2 ≈ fval2     atol=1e-15

    rmin, rmax = big(-1.0), big(1.0)
    imin, imax = big(-1.0), big(1.0)
    g2f = RP.Geometry2Function(rmin, rmax, imin, imax)
    @test eltype(g2f) == BigFloat
end

@testset "RootsAndPoles.jl" begin
    test_GRPFParams()
    test_functionstructs()
end
