distance(p1::Complex, p2::Complex) = hypot(p2 - p1)
distance((p1, p2)::Tuple{Complex, Complex}) = distance(p1, p2)
distance(a::QuadrantPoint, b::QuadrantPoint) = distance(complex(a), complex(b))
distance((a, b)::Tuple{QuadrantPoint, QuadrantPoint}) = distance(complex(a), complex(b))

function sortedge(p1, p2)
    if p1 < p2
        return p1, p2
    else
        return p2, p1
    end
end
sortedge(p::Tuple{Int,Int}) = sortedge(p[1], p[2])

function sorttriangle(u, v, w)
    # Unlike DT.sort_triangle, this does not maintain the orientation of T, e.g.
    # DT.sort_triangle((5, 60, 40)) => (5, 60, 40)
    # DT.sort_triangle((60, 5, 40)) => (5. 40, 60)
    if u < v < w
        return u, v, w
    elseif u < w < v
        return u, w, v
    elseif w < u < v
        return w, u, v
    elseif w < v < u
        return w, v, u
    elseif v < u < w
        return v, u, w
    elseif v < w < u
        return v, w, u
    end
end

"""
    getplotdata(tess, quadrants, phasediffs, g2f)

Return complex coordinates of every tesselation edge and a vector with the corresponding
index specifiying its phase quadrant.

This function is useful for plotting the final tesselation using results from calling
`grpf` with the `PlotData()` argument.

The edge colors are (where Q means "quadrant"):

| Color index |         Meaning         |  Quadrant phase   |
|:-----------:|:-----------------------:|:-----------------:|
|      1      |           Q 1           | 0 ≤ arg f < π/2   |
|      2      |           Q 2           | π/2 ≤ arg f < π   |
|      3      |           Q 3           | π ≤ arg f < 3π/2  |
|      4      |           Q 4           | 3π/2 ≤ arg f < 2π |
|      5      | phase change (boundary) |         -         |
|      6      | candidate edge          |         -         |

# Example

Although not self contained, this example shows how `getplotdata` can be used:

```julia
roots, poles, quads, phasediffs, tess, g2f = grpf(f, origcoords, PlotData());
z, edgecolors = getplotdata(tess, quadrants, phasediffs, g2f);

using Plots
plot(real(z), imag(z), group=edgecolors)
```
"""
function getplotdata(tess, quadrants, phasediffs, g2f)
    edges = [e for e in delaunayedges(tess)]

    edgecolors = fill(5, length(edges))
    for ii in eachindex(edges)
        if phasediffs[ii] == 2  # candidate edges
            edgecolors[ii] = 6
        elseif phasediffs[ii] == 0
            edgecolors[ii] = quadrants[getindex(geta(edges[ii]))]
        end
    end

    edgecolors = repeat(edgecolors, inner=3)  # x, y is xₐ, xᵦ, NaN, repeat

    x, y = getplotxy(edges)
    z = g2f.(x, y)

    I = sortperm(edgecolors)

    return z[I], edgecolors[I]
end

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

        ΔQ = mod(getquadrant(a) - getquadrant(b), 4)  # phase difference
        if ΔQ == 2
            labels[i:i+2] .= 5
        elseif ΔQ == 0
            labels[i:i+2] .= getquadrant(a)
        end
        i += 3
    end
    return xcoords, ycoords, labels
end
