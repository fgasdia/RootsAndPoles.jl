Hfcn(z) = (z + 2)/(z^2 + 1/4)  # zero: -2 and poles: ±im/2
Hfcn_mesh() = rectangulardomain(complex(-3, -1), complex(1, 1), 0.6)

# using GLMakie; Makie.inline!(false)
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

    RP.assignquadrants!(get_points(mesh), Hfcn, false)
    newpts = RP.getquadrant.(RP.each_point(mesh))
    scatter!(ax, reim.(complex.(mesh)), color=colors[newpts], markersize=20)

end

function test_GRPFParams()
    ga = GRPFParams()
    gb = GRPFParams(100, 500000, 3, 1e-9, false)
    gc = GRPFParams(200, 10000, 3, 1e-8, true)

    # equality
    @test ga == gb
    @test ga != gc

    # constructors
    @test GRPFParams(1e-3) == GRPFParams(100, 500000, 3, 1e-3, false)
    @test GRPFParams(1e-3, true) == GRPFParams(100, 500000, 3, 1e-3, true)
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
    RP.assignquadrants!(get_points(tess), Hfcn, true)
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
    RP.assignquadrants!(get_points(tess), Hfcn, false)

    E = Set{Tuple{Int, Int}}()
    RP.candidateedges!(E, tess)

    for e in E
        a, b = RP.get_point(tess, e...)
        @test abs(RP.getquadrant(a) - RP.getquadrant(b)) == 2
    end

    # TODO: test the PlotData argument
end

function test_selectedges!()
    E = Set{Tuple{Int, Int}}()
    selectE = Set{Tuple{Int, Int}}()
    initial_mesh = Hfcn_mesh()
    mesh = RP.QuadrantPoints(RP.QuadrantPoint.(initial_mesh))
    tess = RP.triangulate(mesh)
    RP.assignquadrants!(get_points(tess), Hfcn, false)
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
    tess, E = RP.tesselate!(tess, Hfcn, GRPFParams(10, 5000, 3, 1e-1, false))
    C = RP.contouredges(tess, E)

    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="Re", ylabel="Im")
    # triplot!(ax, tess, triangle_color=colors[RP.getquadrant.(get_points(tess))])

    colors = Makie.colormap("Greens", 2*length(C))
    i = 1
    for c in C
        x1, y1 = reim(complex(get_point(tess, c[1])))
        x2, y2 = reim(complex(get_point(tess, c[2])))
        println(x1," ", y1)
        println(x2," ", y2)
        lines!(ax, [x1, x2], [y1, y2], color=[colors[2i-1], colors[2i]])
        scatter!(ax, [x1, x2], [y1, y2], color=[colors[2i-1], colors[2i]])
        i += 1
    end
    fig
end

function test_evaluateregions()
    mesh = RP.QuadrantPoints(RP.QuadrantPoint.(Hfcn_mesh()))
    tess = RP.triangulate(mesh)
    tess, E = RP.tesselate!(tess, Hfcn, GRPFParams(30, 5000, 3, 1e-3, false))
    C = RP.contouredges(tess, E)

    # XXX: evaluateregions currently fails code_warntype
    regions = RP.evaluateregions!(C, tess)

    fig = Figure()
    ax = Axis(fig[1, 1], xlabel="Re", ylabel="Im")
    # triplot!(ax, tess, triangle_color=colors[RP.getquadrant.(get_points(tess))])
    for r in regions       
        xys = collect(reim.(complex.(get_point(tess, r...))))
        lines!(ax, xys)
        scatter!(ax, xys)
    end
    fig
end

function test_rootsandpoles()
    mesh = RP.QuadrantPoints(RP.QuadrantPoint.(Hfcn_mesh()))
    tess = RP.triangulate(mesh)
    tess, E = RP.tesselate!(tess, Hfcn, GRPFParams(30, 5000, 3, 1e-3, false))
    C = RP.contouredges(tess, E)
    regions = RP.evaluateregions!(C, tess)

    r = regions[1]
    pts = get_point(tess, r...)
    quadrantsequence = [RP.getquadrant(p) for p in pts]

    # Sign flip because `r` are in opposite order of Matlab?
    dq = diff(quadrantsequence)
    for i in eachindex(dq)
        if dq[i] == 3
            dq[i] = -1
        elseif dq[i] == -3
            dq[i] = 1
        elseif abs(dq[i]) == 2
            # ``|ΔQ| = 2`` is ambiguous; cannot tell whether phase increases or
            # decreases by two quadrants
            dq[i] = 0
        end
    end
    q = sum(dq)/4
    z = sum(complex, pts)/length(pts)

    if q > 0
        push!(zroots, z)  # convert in case T isn't Float64
    elseif q < 0
        push!(zpoles, z)
    end
end

@testset "RootsAndPoles.jl" begin
    test_GRPFParams()
    test_assignquadrants!()
    test_candidateedges!()
    test_addzone1node!()
    test_addzone2node!()
end
