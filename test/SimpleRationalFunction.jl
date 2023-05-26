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
    tmap = RP.get_adjacent2vertex(tess)
    for e in unique_pts
        tris = tmap.adjacent2vertex[e]  # ??? What is correct way to handle the map?
        for t in tris
            if t[1] in unique_pts || t[2] in unique_pts
                # t is a zone 2 triangle
                p, q, r = RP.get_point(tess, e, t[1], t[2])
                l1 = RP.distance(p, q)
                l2 = RP.distance(p, r)
                l3 = RP.distance(q, r)
                if max(l1,l2,l3)/min(l1,l2,l3) > parms.skinnytriangle
                    avgnode = (p + q + r)/3
                    # if avgnode !in newnodes
                        println(avgnode)
                        # push!(newnodes, avgnode)
                    # end
                end
            else
                # t is a zone 1 triangle
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
