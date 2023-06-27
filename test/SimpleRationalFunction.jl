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


function edgecolors(tess)
    labels = zeros(Int, 3*num_edges(tess))
    xcoords = Vector{DT.number_type(tess)}(undef, 3*num_edges(tess))
    ycoords = similar(xcoords)

    i = 1
    for e in each_solid_edge(tess)
        a, b = get_point(tess, e[1], e[2])
        xcoords[i] = a.x
        ycoords[i] = a.y
        xcoords[i+1] = b.x
        ycoords[i+1] = b.y
        xcoords[i+2] = NaN
        ycoords[i+2] = NaN

        ΔQ = mod(RP.getquadrant(a) - RP.getquadrant(b), 4)  # phase difference
        if ΔQ == 2
            labels[i:i+2] .= 5
        elseif ΔQ == 0
            labels[i:i+2] .= RP.getquadrant(a)
        end
        i += 3
    end
    return xcoords, ycoords, labels
end


#
mesh = RP.QuadrantPoints(RP.QuadrantPoint.(origcoords))
tess = RP.triangulate(mesh)
RP.assignquadrants!(get_points(tess), simplefcn, false)
pts = RP.getquadrant.(RP.each_point(tess))
ec = edgecolors(tess)

using GLMakie; Makie.inline!(false)
const wc = Makie.wong_colors()
const co = Makie.color
colors = [co("black"), wc[6], wc[7], wc[1], wc[3], wc[4]]
fig = Figure()
ax = Axis(fig, xlabel="Re", ylabel="Im")
limits!(ax, (-1.5, 1.5), (-1.5, 1.5))
lines!(ax, ex, ey, color=colors[ec.+1])

triplot(tess, point_color=colors[pts], show_all_points=true)

# refine!(tess)
tess, E = RP.tesselate!(tess, simplefcn, GRPFParams(1e-6))


C = contouredges(tess, E)
regions = RP.evaluateregions!(C, tess)
zroots, zpoles = RP.rootsandpoles(tess, regions)

@test approxmatch(zroots, matlab_zroots)
@test approxmatch(zpoles, matlab_zpoles)


using Makie.Colors
using GLMakie; Makie.inline!(false)
colors = Makie.wong_colors()[1:4]
fig2 = Figure()
ax = Axis(fig2[1, 1], xlabel="Re", ylabel="Im")
# limits!(ax, (-1.00001,-0.99999), (-0.00001,0.00001))
limits!(ax, (-0.00001, 0.00001), (0.99999, 1.00001))
triplot!(ax, tess, triangle_color=colors[RP.getquadrant.(get_points(tess))])

for r in regions
    cmplxs = [reim(v) for v in complex.(DT.get_point(tess, r...))]

    lines!(ax, cmplxs, color=HSV.(range(0, 360, length(cmplxs)), 50, 50), overdraw=true, linewidth=4)
end
fig2


C = RP.contouredges(tess, E)

# fig = Figure()
# ax = Axis(fig[1, 1], xlabel="Re", ylabel="Im")
# triplot!(ax, tess, triangle_color=colors[RP.getquadrant.(get_points(tess))])


# contouredges
function triangleedges(tess, E)
    # D = Vector{DT.edge_type(tess)}()
    D = Set{DT.edge_type(tess)}()
    for e in E
        # v = get_adjacent(tess, e)  # (e[1], e[2], v) is a positively oriented triangle

        # # All neighboring vertices of the edge
        # v1 = get_neighbours(tess, e[1])
        # v2 = get_neighbours(tess, e[2])
        # # v = symdiff(v1, v2) # if an edge occurs twice...
        # vs = intersect(v1, v2)


        # Get triangles sharing edge `e`
        v1 = get_neighbours(tess, e[1])
        v2 = get_neighbours(tess, e[2])
        vs = intersect(v1, v2)

        for v in vs
            tri = DT.construct_positively_oriented_triangle(tess, e[1], e[2], v)
            i, j, k = indices(tri)
            push!(D, (i, j), (j, k), (k, i))
        end
    end
    # return unique(D)
    return D
end

findall(==(D[60]), D) == [41, 60]
revD = reverse.(D)
findall(==(D[1]), revD) == [4]


function redundant_edge_count(C)
    @assert allunique(C) "Algorithm assumes `C` is unique"
    redundant_count = zeros(Int, length(C))
    for j in eachindex(C)
        if iszero(redundant_count[j])
            # TEMP: replace with more efficient function
            idx = findfirst(==(reverse(C[j])), C)  # there can be no more than 1 match because C isunique
            # length(idx) > 1 && @info idx j
            if isnothing(idx)
                redundant_count[j] = 1
            else
                redundant_count[j] = 2
                redundant_count[idx] = 2
            end
        end
    end

    return redundant_count
end

function isduplicate(C)
    duplicate_edge = Dict{eltype(C),Bool}()
    for c in C
        if reverse(c) in C
            duplicate_edge[c] = true
        else
            duplicate_edge[c] = false
        end
    end

    return duplicate_edge
end
unique_edges(C, duplicate_edge) = Set(c for c in C if !duplicate_edge[c])

# unique_edges(C, redundant_count) = [c for (i, c) in pairs(C) if redundant_count[i] == 1]

function contouredges(tess, E)
    D = triangleedges(tess, E)
    duplicates = isduplicate(D)
    C = unique_edges(D, duplicates)
    return C
end


C = contouredges(tess, E)
regions = RP.evaluateregions!(C, tess)
zroots, zpoles = RP.rootsandpoles(tess, regions)


    for c in C
        x1, y1 = reim(complex(DT.get_point(tess, c[1])))
        x2, y2 = reim(complex(DT.get_point(tess, c[2])))
        lines!(ax, [x1, x2], [y1, y2], color="red", overdraw=true, linewidth=4)
        # scatter!(ax, [x1, x2], [y1, y2], color=[colors[2i-1], colors[2i]])
        sleep(0.5)
        fig2
    end
    fig2
    




# evaluateregions
c = first(C)
regions = [[c[1]]]
refidx = c[2]
numregions = 1
delete!(C, c)

nextedges = Vector{DT.edge_type(tess)}()
while length(C) > 0
    # Writing out the for loop avoids a Core.box closure issue with `refidx`
    # nextedges = [c for c in C if refidx in c]
    for c in C
        if c[1] == refidx
            push!(nextedges, c)
        end
    end

    if isempty(nextedges)
        push!(regions[numregions], refidx)

        c = first(C)
        numregions += 1
        push!(regions, [c[1]])
        refidx = c[2]
        delete!(C, c)
    else
        if length(nextedges) > 1
            previdx = last(regions[numregions])
            c = findnextedge(tess, previdx, refidx, nextedges)
        else
            c = only(nextedges)
        end
        push!(regions[numregions], c[1])
        refidx = c[2]
        delete!(C, c)
    end
    empty!(nextedges)
end
push!(regions[numregions], refidx)














for r in regions
    cmplxs = [reim(v) for v in complex.(DT.get_point(tess, r...))]

    lines!(ax, cmplxs, color=HSV.(range(0, 360, length(cmplxs)), 50, 50), overdraw=true, linewidth=4)
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
end

fig2







zroots = Vector{complex(DT.number_type(tess))}()
zpoles = similar(zroots)
for r in regions
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
