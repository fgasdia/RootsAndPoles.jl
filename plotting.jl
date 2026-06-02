using GLMakie
GLMakie.activate!()
import GLMakie.ColorTypes.ColorTypes as CT

using DelaunayTriangulation
using RootsAndPoles
const RP = RootsAndPoles

function formatquadrants(mesh)
    xs = Vector{Float64}()
    ys = similar(xs)
    cs = Vector{Int}()
    for e in RP.each_solid_edge(mesh)
        u, v = edge_vertices(e)
        x, y = RP.coord(mesh, u)
        p, q = RP.coord(mesh, v)

        qu, qv = RP.quadrant(mesh, u, v)
        ΔQ = abs(qu - qv)
        ΔQ == 2 ? c = 5 :
            ΔQ == 0 ? c = RP.quadrant(mesh, u) :
                c = 0

        push!(xs, x, p)
        push!(ys, y, q)
        push!(cs, c, c)
    end
    return xs, ys, cs
end

function formatcontours(mesh, contours)
    xs = Vector{Float64}()
    ys = similar(xs)
    for c in contours
        for e in c
            u, v = edge_vertices(e)
            x, y = RP.coord(mesh, u)
            p, q = RP.coord(mesh, v)

            push!(xs, x, p)
            push!(ys, y, q)
        end
    end
    return xs, ys
end

function formatcom(mesh, contours)
    xs = Vector{Float64}()
    ys = similar(xs)
    for c in contours
        z = RP.centerofmass(mesh, c)
        x, y = reim(z)
        push!(xs, x)
        push!(ys, y)
    end
    return xs, ys
end

function plotquadrants(mesh; title)
    fig = Figure(size=(650, 650))
    wh = (width=580, height=580)
    ax = Axis(fig[1,1]; wh..., title, xlabel="Re(z)", ylabel="Im(z)")
    xs, ys, cs = formatquadrants(mesh)
    tab10cmap = to_colormap(:tab10)
    quadrantcolormap = [CT.RGBA(0,0,0,1), tab10cmap[4], tab10cmap[2], tab10cmap[3], tab10cmap[1], tab10cmap[5]]
    lwidths = ones(length(cs))
    candlw = 3
    lwidths[cs .== 5] .= candlw
    ls = linesegments!(ax, xs, ys, color=cs, colormap=quadrantcolormap, colorrange=(0,5),
                       linewidth = lwidths,
                       label = [label => (; color=i, linewidth=i < 5 ? 1 : candlw)
                                for (i, label) in enumerate(["0 ≤ arg f < π/2",
                                "π/2 ≤ arg f < π", "π ≤ arg f < 3π/2", "3π/2 ≤ arg f < 2π",
                                "candidate"])])
    axislegend(ax, position=:rc)
    display(GLMakie.Screen(), fig)

    return fig, ax
end

function plotrefinement(mesh, type)
    xs = Vector{Float64}()
    ys = similar(xs)
    cs = Vector{Int}()
    for e in each_solid_edge(mesh)
        u, v = edge_vertices(e)
        x, y = get_point(mesh, u)
        p, q = get_point(mesh, v)

        if e in keys(type)
            c = Int(type[e])
        elseif reverse(e) in keys(type)
            c = Int(type[reverse(e)])
        else
            c = 0
        end

        push!(xs, x, p)
        push!(ys, y, q)
        push!(cs, c, c)
    end
    return xs, ys, cs
end

function plotquadrants(iterations::MeshIterations, i, plotcontours=true)
    mesh = iterations.mesh[i]
    fig, ax = plotquadrants(mesh; title="iter $i - $(iterations.refinement_mode[i])")
    if plotcontours && any(==(RP.SpecialEdge.Candidate), values(iterations.splitedges[i]))
        contours = RP.buildregions(mesh)
        xs, ys = formatcontours(mesh, contours)
        linesegments!(ax, xs, ys, color=:brown, linewidth=3)

        xs, ys = formatcom(mesh, contours)
        scatter!(ax, xs, ys, color=:black, markersize=10)
    end
    return fig, ax
end
