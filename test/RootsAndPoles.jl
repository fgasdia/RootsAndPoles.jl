#==
o Structs
    x ComplexMesh
    x FinderParams
    o SpecialEdge ?
    o MeshIterations: are lengths correct?
    o PreviousIteration ?
    x GradientAnalysis : empty!

o interfaces
    o rootsandpoles
        o ReturnMultiplicity
        o SkinnyMode
        o MeshIterations
        o test # of maxnodes is obeyed, maxiter, etc

o functions
    o trianglegradient: confirm grad algorithm to matlab expression?
    o analyzegradients!: make a plot to confirm it's working?
    o evalfcn!: confirm that ComplexMesh is updating, that threads are used, etc


# tests
o contouredges, increment_solid!, and edgekey: time this approach vs minmax on edge
==#


# TODO: Does max length keep getting shorter over iterations?
# TODO: large and small tolerance
# TODO: Write test comparing skinnymode criteria
# TODO: Test with dense initial coords and sparse initial coords
# TODO: Make sure it works to take the modified mesh and use it as starting point for new rootsandpoles? need to set startidx

@testset "FinderParams" begin
    @test isconcretetype(FinderParams)

    pa = FinderParams()
    pb = FinderParams(1e-9, 1000, 2000, 10000, 1)
    pc = FinderParams(tol=1e-8)

    @test pa == pb
    @test pa != pc
end

@testset "ComplexMesh" begin
    coords = [-1-1im, 1-1im, 1+1im, -1+1im]
    scoords = reim.(float(coords))
    a = @inferred ComplexMesh(coords)
    b = @inferred ComplexMesh(coords; rng=Random.default_rng())
    c = @inferred ComplexMesh(scoords)

    @test a.coords == b.coords
    @test a.coords == c.coords

    @test num_solid_vertices(a) == DT.num_solid_vertices(a.tri)
    @test length(coords) == length(a.coords)
    @test length(a.coords) == DT.num_solid_vertices(a.tri)
    @test length(a.coords) == length(a.fvals)

    @test lastindex(a) == lastindex(a.coords)
    
    d = copy(a)
    @test d.coords == a.coords
    @test d.fvals == a.fvals
    @test d.tri == a.tri
    @test d.rng == a.rng

    push!(d.coords, (0.0, 0.0))
    @test a.coords != d.coords

    push!(c.coords, (0.0, 0.0))
    @test c.coords == scoords  # modifies scoords in-place

    # types
    @test DT.triangle_type(a.tri) == RP.TRIANGLETYPE
    @test DT.edge_type(a.tri) == RP.EDGETYPE
    @test DT.number_type(a.tri) == Float64

    # coord
    c1 = @inferred RP.coord(a, 1)
    c2 = @inferred RP.coord(a, 1, 2)
    coordtype = NTuple{2, Float64}
    @test c1 isa coordtype
    @test c2 isa Tuple{coordtype, coordtype}

    for i in eachindex(coords)
        cai = RP.coord(a,i)
        @test cai == reim(float(coords[i]))
        @test cai == DT.get_point(a.tri, i)
    end

    # fval
    f1 = @inferred RP.fval(a, 1)
    f2 = @inferred RP.fval(a, 1, 2)
    ftype = ComplexF64
    @test f1 isa ftype
    @test f2 isa Tuple{ftype, ftype}
    fs = @inferred RP.fvals(a)
    @test fs == a.fvals
    RP.fval(a, 1) == 0.123 ? (testval = 0.222) : (testval = 0.123)
    @inferred RP.fval!(a, 1, testval)
    @test RP.fval(a, 1) == testval

    # quadrant
    q1 = @inferred RP.quadrant(b, 1)
    q2 = @inferred RP.quadrant(b, 1, 2)
    @test length(q1) == 1
    @test length(q2) == 2
    @test q1 == RP.quadrant(RP.fval(b, 1))

    # edge_length
    firstedge = first(each_solid_edge(a.tri))
    e1 = @inferred RP.edge_length(a, firstedge)
    @test e1 == DT.edge_length(a.tri, firstedge)
    @test e1 isa Float64

    # each_solid_...
    @inferred RP.each_solid_edge(a)
    @inferred RP.each_solid_triangle(a)
    
    # get_...
    @inferred RP.get_adjacent(a, firstedge)
    @inferred RP.get_adjacent(a, firstedge[1], firstedge[2])
    @inferred RP.get_neighbours(a, 1)

    # midpoint
    m = @inferred RP.midpoint(a, 1, 2)
    @test m == (0.0, -1.0)  # between (-1, -1) and (1, -1)

    # addpoints
    @test length(scoords) > num_solid_vertices(c)  # we pushed to scoords above
    @inferred RP.addpoints!(c)
    for i in eachindex(scoords)
        @test DT.get_point(c.tri, i) == scoords[i]
    end

    # rng
    rng1 = Random.MersenneTwister()
    rng2 = Random.Xoshiro()
    r1 = ComplexMesh(coords; rng=rng1)
    r2 = ComplexMesh(coords; rng=rng2)
    @test r1.rng == rng1
    @test r2.rng == rng2

    @testset "evalfcn!" begin
        fun = abs2
        @inferred RP.evalfcn!(fun, a, 1)

        @test RP.fval(a, 2) == fun(coords[2])
        @test RP.fval(a, 3) == fun(coords[3])

        @test RP.quadrant(a, 2) == RP.quadrant(fun(coords[2]))
        @test RP.quadrant(a, 3) == RP.quadrant(fun(coords[3]))

        fun2(z) = 3
        RP.evalfcn!(fun2, a, 3)  # only apply fun2 to indices 3 and 4

        @test RP.fval(a, 2) == fun(complex(coords[2]))
        @test RP.fval(a, 3) == fun2(complex(coords[3]))
    end
end

@testset "GradientAnalysis" begin
    @test isconcretetype(RP.GradientAnalysis)

    ga = RP.GradientAnalysis(
        Dict{RP.TRIANGLETYPE,SVector{2,Float64}}(),
        Vector{Float64}(),
        Vector{Float64}()
    )
    push!(ga.gradients, 
        (1,2,3)=>SVector(0.0,0.0), (4,5,6)=>SVector(0.0, 1.0), (7,8,9)=>SVector(1.0,1.5)
    )
    push!(ga.indicators, 0.0, 0.1, 0.2, 0.3)
    push!(ga.remainingedgelengths, 1.0, 2.0, 3.0, 4.0)

    @test length(ga.indicators) == 4
    @test length(ga.remainingedgelengths) == 4
    @test length(ga.gradients) == 3

    @inferred empty!(ga)
    for f in fieldnames(typeof(ga))
        @test isempty(getfield(ga, f))
    end
end

########

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

@testset "deduplication" begin
    @test RP.toldigits(1e-9) == 8
    @test RP.toldigits(5e-9) == 8
    @test RP.toldigits(1e-3) == 2
    @test RP.toldigits(1) == -1
    @test RP.toldigits(10) == -2

    v = [1.332, 1.331, 1.395, 1.320, 1.301]
    v2 = copy(v)
    RP.dedupe!(v; digits=RP.toldigits(1e-3))
    RP.dedupe!(v2; tol=1e-3)
    @test v == [1.332, 1.395, 1.320, 1.301]
    @test v == v2
end

@testset "RootsAndPoles.jl" begin
   
    
end
