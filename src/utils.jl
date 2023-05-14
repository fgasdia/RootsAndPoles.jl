distance(p1::Complex, p2::Complex) = hypot(p2 - p1)
distance((p1, p2)::Tuple{Complex, Complex}) = distance(p1, p2)
distance(a::QuadrantPoint, b::QuadrantPoint) = distance(complex(a), complex(b))

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
