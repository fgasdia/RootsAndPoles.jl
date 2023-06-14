Hfcn(z) = (z + 2)/(z^2 + 1/4)  # zero: -2 and poles: ±im/2
Hfcn_mesh() = rectangulardomain(complex(-3, -1), complex(1, 1), 0.6)

using GLMakie; Makie.inline!(false)
function Hfcn_plot()
    initial_mesh = Hfcn_mesh()
    mesh = RP.QuadrantPoints(RP.QuadrantPoint.(initial_mesh))
    tess = triangulate(mesh)
    RP.assignquadrants!(get_points(tess), Hfcn, false)
    pts = RP.getquadrant.(RP.each_point(mesh))
    colors = Makie.wong_colors()[1:4]

    
    limits!(ax, -3.5, 1.5, -1.2, 1.2)
    scatter!(ax, reim.(initial_mesh), color=colors[pts], markersize=10)
    elements = [PolyElement(polycolor = colors[i]) for i in 1:4]
    Legend(f[1,2], elements, string.(1:4), label="Quadrant")
    f

    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="Re", ylabel="Im")
    triplot!(ax, tess, triangle_color=colors[pts], point_color=colors[pts])


    E = Set{Tuple{Int, Int}}()
    selectE = Set{Tuple{Int, Int}}()
    RP.assignquadrants!(get_points(tess), Hfcn, false)
    RP.candidateedges!(E, tess)
    RP.selectedges!(selectE, tess, E, GRPFParams().tolerance)
    unique_pts = Set(Iterators.flatten(selectE))
    RP.splittriangles!(tess, unique_pts, GRPFParams().tolerance, GRPFParams().skinnytriangle)

    RP.assignquadrants!(get_points(tess), Hfcn, false)
    triplot!(ax, tess, triangle_color=colors[RP.getquadrant.(get_points(tess))], point_color=colors[pts])



    scatter!(ax, reim.(initial_mesh), color=colors[pts], markersize=10)
    elements = [PolyElement(polycolor = colors[i]) for i in 1:4]
    Legend(f[1,2], elements, string.(1:4), label="Quadrant")
    f

    RP.assignquadrants!(mesh, Hfcn, false)
    newpts = RP.getquadrant.(RP.each_point(mesh))
    scatter!(ax, reim.(complex.(mesh)), color=colors[newpts], markersize=20)

end

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
    pts = Hfcn_mesh()

    # multithreading = false
    mesh = RP.QuadrantPoints(RP.QuadrantPoint.(pts))
    tess = triangulate(mesh)
    @test iszero(RP.getquadrant.(get_points(tess)))  
    RP.assignquadrants!(get_points(tess), Hfcn, false)
    @test RP.getquadrant.(mesh.points) == RP.quadrant.(Hfcn.(pts))

    # multithreading = true
    mesh = RP.QuadrantPoints(RP.QuadrantPoint.(pts))
    tess = triangulate(mesh)
    RP.assignquadrants!(tess, Hfcn, true)
    @test RP.getquadrant.(get_points(tess)) == RP.quadrant.(Hfcn.(pts))
end

# function test_newset()
#     S = Set{Tuple{Int,Int}}()
#     for _ in 1:100
#         s = (rand(1:100), rand(1:100))
#         push!(S, s)
#     end
#     return S
# end

# S = Set{Tuple{Int,Int}}(tuple.(rand(1:100, 100), rand(1:100, 100)))
# function test_emptyset(S)
#     empty!(S)
#     for _ in 1:100
#         s = (rand(1:100), rand(1:100))
#         push!(S, s)
#     end
#     return S
# end

function test_candidateedges!()  
    initial_mesh = Hfcn_mesh()
    mesh = RP.QuadrantPoints(RP.QuadrantPoint.(initial_mesh))
    tess = RP.triangulate(mesh)
    RP.assignquadrants!(tess, Hfcn, false)

    E = Set{Tuple{Int, Int}}()
    RP.candidateedges!(E, tess)

    for e in E
        a, b = RP.get_point(tess, e...)
        @test abs(RP.getquadrant(a) - RP.getquadrant(b)) == 2
    end

    # TODO: test the PlotData argument

    unique_pts = Set(Iterators.flatten(selectE))
    empty!(mesh_points)  # XXX XXX Doing this clears the points in tess!! Look at DT animation for adding
    RP.splittriangles!(mesh_points, tess, unique_pts, GRPFParams().tolerance, GRPFParams().skinnytriangle)
end

function test_selectedges!()
    E = Set{Tuple{Int, Int}}()
    selectE = Set{Tuple{Int, Int}}()
    initial_mesh = Hfcn_mesh()
    mesh = RP.QuadrantPoints(RP.QuadrantPoint.(initial_mesh))
    tess = RP.triangulate(mesh)
    RP.assignquadrants!(tess, Hfcn, false)
    RP.candidateedges!(E, tess)

    RP.selectedges!(selectE, tess, E, GRPFParams().tolerance)
    @test selectE == E

    RP.selectedges!(selectE, tess, E, 1)
    @test isempty(selectE)
end

function test_addzone1node!()
    initial_mesh = [0+0im, 1+0im, 1+1im, 0+1im]
    mesh = RP.QuadrantPoints(RP.QuadrantPoint.(initial_mesh))
    tess = RP.triangulate(mesh)
    RP.addzone1node!(tess, 0+0im, 1+1im, GRPFParams().tolerance)
    @test 0.5+0.5im in complex.(get_points(tess))
end

function test_addzone2node!()
    initial_mesh = [0+0im, 1+0im, 1+1im, 0+1im]
    mesh = RP.QuadrantPoints(RP.QuadrantPoint.(initial_mesh))
    tess = RP.triangulate(mesh)
    RP.addzone2node!(tess, 0+0im, 1+0im, 1+1im, GRPFParams().tolerance)
    @test (2/3)+(1/3)im in complex.(get_points(tess))
end

function test_tesselate!()
    mesh = RP.QuadrantPoints(RP.QuadrantPoint.(Hfcn_mesh()))
    tess = RP.triangulate(mesh)
    params = GRPFParams(5, 5000, 3, 1e-9, false)
    tess, E = RP.tesselate!(tess, Hfcn, params)
    
    colors = Makie.wong_colors()[1:4]
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="Re", ylabel="Im")
    triplot!(ax, tess, triangle_color=colors[RP.getquadrant.(get_points(tess))])
end

function test_contouredges()
    mesh = RP.QuadrantPoints(RP.QuadrantPoint.(Hfcn_mesh()))
    tess = RP.triangulate(mesh)
    tess, E = RP.tesselate!(tess, Hfcn, GRPFParams(30, 5000, 3, 1e-3, false))
    C = RP.contouredges(tess, E)

    colors = Makie.wong_colors()[1:4]
    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="Re", ylabel="Im")
    # triplot!(ax, tess, triangle_color=colors[RP.getquadrant.(get_points(tess))])
    for c in C
        x1, y1 = reim(complex(get_point(tess, c[1])))
        x2, y2 = reim(complex(get_point(tess, c[2])))
        println(x1," ", y1)
        println(x2," ", y2)
        lines!(ax, [x1, x2], [y1, y2])
    end
end

function test_evaluateregions()
    mesh = RP.QuadrantPoints(RP.QuadrantPoint.(Hfcn_mesh()))
    tess = RP.triangulate(mesh)
    tess, E = RP.tesselate!(tess, Hfcn, GRPFParams(30, 5000, 3, 1e-3, false))
    C = RP.contouredges(tess, E)

    #
    numregions = 1

    c = first(C)
    regions = [[c[1]]]
    refpt = c[2]
    delete!(C, c)

    nextedges = collect(Iterators.filter(v->v[1] == refpt, C))
    prevpt = regions[numregions][end]
    RP.findnextpt(tess, prevpt, refpt, nextedges)
    #

    regions = RP.evaluateregions!(C, tess)
end

@testset "RootsAndPoles.jl" begin
    test_GRPFParams()
    test_functionstructs()
    test_assignquadrants!()
    test_candidateedges!()
    test_addzone1node!()
    test_addzone2node!()
end
