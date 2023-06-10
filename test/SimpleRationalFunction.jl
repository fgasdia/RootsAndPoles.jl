function simplefcn(z)
    w = (z - 1)*(z - im)^2*(z + 1)^3/(z + im)
end

# Analysis parameters
xb = -2  # real part begin
xe = 2  # real part end
yb = -2  # imag part begin
ye = 2  # imag part end
r = 0.1  # initial mesh step



###
# TEMP
const RP = RootsAndPoles
const DT = RP.DelaunayTriangulation

parms = GRPFParams()
initial_mesh = RP.rectangulardomain(complex(xb, yb), complex(xe, ye), r)

mesh_points = RP.QuadrantPoints(RP.QuadrantPoint.(initial_mesh))
RP.assignquadrants!(mesh_points, simplefcn, false)

E = Set{Tuple{Int, Int}}()  # edges
selectE = Set{Tuple{Int, Int}}()
tess = DT.triangulate(mesh_points)
RP.candidateedges!(E, tess)



H(z) = (z + 2)/(z^2 + 1/4)  # zero: -2 and poles: ±im/2
initial_mesh = rectangulardomain(complex(-3, -1), complex(1, 1), 0.6)

using GLMakie
scatter(reim.(initial_mesh))

mesh = RP.QuadrantPoints(RP.QuadrantPoint.(initial_mesh))
RP.assignquadrants!(mesh, H, false)
pts = RP.getquadrant.(RP.each_point(mesh))
colors = Makie.wong_colors()[1:4]
f = Figure()
ax = Axis(f[1, 1], xlabel="Re", ylabel="Im")
limits!(ax, -3.5, 1.5, -1.2, 1.2)
scatter!(ax, reim.(initial_mesh), color=colors[pts], markersize=10, label="1")
elements = [PolyElement(polycolor = colors[i]) for i in 1:4]
Legend(f[1,2], elements, string.(1:4))
f

RP.assignquadrants!(mesh_points, H, false)
newpts = RP.getquadrant.(RP.each_point(mesh_points))
scatter!(ax, reim.(complex.(mesh_points)), color=colors[newpts], markersize=20, label="2")


maxElength = 0.0
for e in E
    d = RP.distance(RP.get_point(tess, e...))
    println(d)
    if d > parms.tolerance
        push!(selectE, e)
        if d > maxElength
            maxElength = d
        end
    end
end

unique_pts = Set(Iterators.flatten(selectE))

function zone2a(tess, unique_pts, parms)
    for T in RP.each_triangle(tess)
        i, j, k = RP.indices(T)
        zone = 0
        if i in unique_pts
            zone += 1
        end
        if j in unique_pts
            zone += 1
        end
        if k in unique_pts
            zone += 1
        end
    
        if zone > 1
            p, q, r = RP.get_point(tess, i, j, k)
            l1 = RP.distance(p, q)
            l2 = RP.distance(p, r)
            l3 = RP.distance(q, r)
            if max(l1,l2,l3)/min(l1,l2,l3) > parms.skinnytriangle
                avgnode = (p + q + r)/3
                println(avgnode)
                # push!(newnodes, avgnode)
            end
        end
    end 
end    
@benchmark zone2a(tess, unique_pts, parms)

# zone2b appears much faster than zone2a
function zone2b(tess, unique_pts, parms)
    triangles = Set{DT.triangle_type(tess)}()
    edges = Set{DT.edge_type(tess)}()
    for p1 in unique_pts
        adj2v = RP.get_adjacent2vertex(tess, p1) 
        for (p2, p3) in adj2v
            sortedtri = DT.sort_triangle((p1, p2, p3))
            if sortedtri in triangles
                continue  # this triangle has already been visited
            end
            push!(triangles, sortedtri)
            p, q, r = RP.get_point(tess, p1, p2, p3)

            if p2 in unique_pts || p3 in unique_pts
                # (p1, p2, p3) is a zone 1 triangle
                # Add a new node at the midpoint of each edge of (p1, p2, p3)
                if !(sort_edge(p1, p2) in edges)
                    addzone1node!(newnodes, p, q, parms.tolerance)
                    push!(edges, sort_edge(p1, p2))
                end
                if !(sort_edge(p1, p3) in edges)
                    addzone1node!(newnodes, p, r, parms.tolerance)
                    push!(edges, sort_edge(p1, p3))
                end
                if !(sort_edge(p2, p3) in edges)
                    addzone1node!(newnodes, q, r, parms.tolerance)
                    push!(edges, sort_edge(p2, p3))
                end
            else
                # (p1, p2, p3) is a zone 2 triangle
                # Add a new node at the average of (p1, p2, p3) 
                addzone2node!(newnodes, p, q, r, parms.skinnytriangle)
            end
        end
    end
end
@benchmark zone2b(tess, unique_pts, parms)

###


origcoords = rectangulardomain(complex(xb, yb), complex(xe, ye), r)

# matlab results from https://github.com/PioKow/GRPF for comparison
matlab_zroots = [-0.999999999951224 - 0.000000000028656im,
                  0.000000000253637 + 1.000000000074506im,
                  1.000000000317046 - 0.000000000062088im]

matlab_zpoles = [0.000000000380455 - 0.999999999701977im]

zroots, zpoles = grpf(simplefcn, origcoords)

@test length(zroots) == 3
@test length(zpoles) == 1

@test approxmatch(zroots, matlab_zroots)
@test approxmatch(zpoles, matlab_zpoles)

pzroots, pzpoles, quadrants, phasediffs, tess, g2f = grpf(simplefcn, origcoords, PlotData())

@test approxmatch(pzroots, matlab_zroots)
@test approxmatch(pzpoles, matlab_zpoles)

# Test with big origcoords
xb = big"-2"  # real part begin
xe = big"2"  # real part end
yb = big"-2"  # imag part begin
ye = big"2"  # imag part end
r = big"0.1"  # initial mesh step

origcoords = rectangulardomain(complex(xb, yb), complex(xe, ye), r)

bzroots, bzpoles = grpf(simplefcn, origcoords)

@test all(isa.(bzroots, Complex{BigFloat}))
@test all(isa.(bzpoles, Complex{BigFloat}))

@test approxmatch(bzroots, zroots)
@test approxmatch(bzpoles, zpoles)
