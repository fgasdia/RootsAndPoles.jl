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

function test_assignquadrants!()
    pts = [complex(1.5, 0.2), complex(-0.4, 2.0), complex(0.9, -0.3), complex(-3.4, -2.0)]
    
    # multithreading = false
    mesh = RP.QuadrantPoints(RP.QuadrantPoint.(pts))
    @test iszero(RP.getquadrant.(mesh.points))
    RP.assignquadrants!(mesh, exp, false)
    @test RP.getquadrant.(mesh.points) == RP.quadrant.(exp.(pts))

    # multithreading = true
    mesh = RP.QuadrantPoints(RP.QuadrantPoint.(pts))
    RP.assignquadrants!(mesh, exp, true)
    @test RP.getquadrant.(mesh.points) == RP.quadrant.(exp.(pts))
end

# function test_newset()
#     S = Set{Tuple{Int,Int}}()
#     for _ in 1:10_000
#         s = (rand(1:1000), rand(1:1000))
#         push!(S, s)
#     end
#     return S
# end

# S = Set{Tuple{Int,Int}}(tuple.(rand(1:1000, 10000), rand(1:1000, 10000)))
# function test_emptyset(S)
#     empty!(S)
#     for _ in 1:10_000
#         s = (rand(1:1000), rand(1:1000))
#         push!(S, s)
#     end
#     return S
# end

function test_candidateedges!()
    H(z) = (z + 2)/(z^2 + 1/4)  # zero: -2 and poles: ±im/2
    
    initial_mesh = rectangulardomain(complex(-3, -1), complex(1, 1), 0.6)
    mesh_points = RP.QuadrantPoints(RP.QuadrantPoint.(initial_mesh))
    RP.assignquadrants!(mesh_points, H, false)

    E = Set{Tuple{Int, Int}}()
    tess = RP.triangulate(mesh_points)
    RP.candidateedges!(E, tess)

    for e in E
        a, b = RP.get_point(tess, e...)
        @test abs(RP.getquadrant(a) - RP.getquadrant(b)) == 2
    end

    # selectE = Set{Tuple{Int, Int}}()
    # RP.selectedges!(selectE, tess, E, GRPFParams().tolerance)

    unique_pts = Set(Iterators.flatten(selectE))
    empty!(mesh_points)  # XXX XXX Doing this clears the points in tess!! Look at DT animation for adding
    RP.splittriangles!(mesh_points, tess, unique_pts, GRPFParams().tolerance, GRPFParams().skinnytriangle)
end

@testset "RootsAndPoles.jl" begin
    test_GRPFParams()
    test_functionstructs()
    test_assignquadrants!()
    test_candidateedges!()
end
