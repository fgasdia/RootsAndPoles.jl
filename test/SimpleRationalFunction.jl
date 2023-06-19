function simplefcn(z)
    w = (z - 1)*(z - im)^2*(z + 1)^3/(z + im)
end

# Analysis parameters
xb = -2  # real part begin
xe = 2  # real part end
yb = -2  # imag part begin
ye = 2  # imag part end
r = 0.1  # initial mesh step

origcoords = rectangulardomain(complex(xb, yb), complex(xe, ye), r)

# matlab results from https://github.com/PioKow/GRPF for comparison
matlab_zroots = [-0.999999999951224 - 0.000000000028656im,
                  0.000000000253637 + 1.000000000074506im,
                  1.000000000317046 - 0.000000000062088im]

matlab_zpoles = [0.000000000380455 - 0.999999999701977im]


#
mesh = RP.QuadrantPoints(RP.QuadrantPoint.(origcoords))
tess = RP.triangulate(mesh)
tess, E = RP.tesselate!(tess, simplefcn, GRPFParams(1e-4))

colors = Makie.wong_colors()[1:4]
fig = Figure()
ax = Axis(fig[1, 1], xlabel="Re", ylabel="Im")
triplot!(ax, tess, triangle_color=colors[RP.getquadrant.(get_points(tess))])

C = RP.contouredges(tess, E)

# fig = Figure()
# ax = Axis(fig[1, 1], xlabel="Re", ylabel="Im")
# triplot!(ax, tess, triangle_color=colors[RP.getquadrant.(get_points(tess))])


# contouredges

C = Set{DT.edge_type(tess)}()
# D = Vector{DT.edge_type(tess)}()
for e in E
    # v = get_adjacent(tess, e)  # (e[1], e[2], v) is a positively oriented triangle

    # All neighboring vertices of the edge
    v1 = get_neighbours(tess, e[1])
    v2 = get_neighbours(tess, e[2])
    # v = symdiff(v1, v2) # if an edge occurs twice...
    vs = intersect(v1, v2)

    # v1 = get_adjacent2vertex(tess, e[1])
    # v2 = get_adjacent2vertex(tess, e[2])

    for v in vs
        if (e[1], v) in C
            delete!(C, (e[1], v))
        elseif (v, e[1]) in C
            delete!(C, (v, e[1]))
        else
            push!(C, (e[1], v))
        end

        if (e[2], v) in C
            delete!(C, (e[2], v))
        elseif (v, e[2]) in C
            delete!(C, (v, e[2]))
        else
            push!(C, (e[2], v))
        end
    end
end

    tri = triangulate(RP.QuadrantPoints(collect(get_point(tess, v...))))
    D = get_convex_hull(tri)
    triplot(tri, convex_hull_linewidth=6)

    hull = convex_hull(RP.QuadrantPoints(collect(get_point(tess, v...))))
    # DT.num_points(t::NTuple) = length(t)
    # convex_hull(reim.(collect(complex.(get_point(tess, v...)))))


    # If an edge occurs twice, that is because it is an edge shared by multiple triangles
    # and by definition is not an edge on the boundary of the candidate region.
    # Therefore, if an edge we come to is already in C, delete it.
    # if (v, e[1]) in D
    #     # delete!(C, (v, e[1]))
    # elseif (e[1], v) in D
    #     # We need to check both edge orientations (a, b) and (b, a)
    #     # delete!(C, (e[1], v))
    # else
        push!(D, (v, e[1]))
    # end

    # if (v, e[2]) in D
    #     # delete!(C, (v, e[2]))
    # elseif (e[2], v) in D
    #     # delete!(C, (e[2], v))
    # else
        push!(D, (e[2], v))
    # end

    # if e in D
    #     # delete!(C, e)
    # elseif (e[2], e[1]) in D
    #     # delete!(C, (e[2], e[1]))
    # else
        push!(D, e)
    # end
end


for d in D
    x1, y1 = reim(complex(DT.get_point(tess, d[1])))
    x2, y2 = reim(complex(DT.get_point(tess, d[2])))
    # println(x1," ", y1)
    # println(x2," ", y2)
    lines!(ax, [x1, x2], [y1, y2], color="lime", overdraw=true, linewidth=4)
    # scatter!(ax, [x1, x2], [y1, y2], color=[colors[2i-1], colors[2i]])
    i += 1
end
#

colors = Makie.colormap("Reds", 2*length(C))
i = 1
for c in C
    x1, y1 = reim(complex(DT.get_point(tess, c[1])))
    x2, y2 = reim(complex(DT.get_point(tess, c[2])))
    # println(x1," ", y1)
    # println(x2," ", y2)
    lines!(ax, [x1, x2], [y1, y2], color="red", overdraw=true, linewidth=4)
    # scatter!(ax, [x1, x2], [y1, y2], color=[colors[2i-1], colors[2i]])
    i += 1
end


# scatter!(ax, reim(collect(complex.(DT.get_point(tess, Iterators.flatten(E)...))))...,
    # color="white", markersize=10,overdraw=true)

for e in E
    x1, y1 = reim(complex(DT.get_point(tess, e[1])))
    x2, y2 = reim(complex(DT.get_point(tess, e[2])))
    # println(x1," ", y1)
    # println(x2," ", y2)
    lines!(ax, [x1, x2], [y1, y2], color="white", overdraw=true, linewidth=4)
    scatter!(ax, [x1, x2], [y1, y2], color="white", overdraw=true, markersize=10)
    i += 1
end

fig
#




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
