# GRPF.jl: Global complex Roots and Poles Finding in Julia

[![Build Status](https://travis-ci.com/fgasdia/GRPF.jl.svg?branch=master)](https://travis-ci.com/fgasdia/GRPF.jl) [![Build status](https://ci.appveyor.com/api/projects/status/megpgn8l1ej5m3ww?svg=true)](https://ci.appveyor.com/project/fgasdia/grpf-jl) [![DOI](https://zenodo.org/badge/154031378.svg)](https://zenodo.org/badge/latestdoi/154031378)

A Julia implementation of [GRPF](https://github.com/PioKow/GRPF) by Piotr Kowalczyk.

## Description

GRPF attempts to **find all the zeros and poles of a complex valued function
with complex arguments in a fixed region**. These types of problems are frequently encountered in electromagnetics, but the algorithm can also be used for similar problems in e.g. optics, acoustics, etc.

GRPF first samples the function on a triangular mesh through Delaunay triangulation. Candidate regions to search for roots and poles are determined and the discretized [Cauchy's argument principle](https://en.wikipedia.org/wiki/Argument_principle) is applied _without needing the derivative of the function or integration over the contour_. To improve the accuracy of the results, a self-adaptive mesh refinement occurs inside the identified candidate regions.

![simplefcn](simplefcn.svg)

## Usage

### Installation

```julia
]add https://github.com/fgasdia/GRPF.jl
```

### Example Problem

Consider a simple rational (complex) argument function `simplefcn(z)` for which we seek roots and poles.
```julia
function simplefcn(z)
    w = (z - 1)*(z - im)^2*(z + 1)^3/(z + im)
end
```

Next, we define parameters for the initial complex domain over which we would like to search.
```julia
xb = -2  # real part begin
xe = 2  # real part end
yb = -2  # imag part begin
ye = 2  # imag part end
r = 0.1  # initial mesh step
tolerance = 1e-9
```

This package includes convenience functions for rectangular and disk shaped domains, but any "shape" can be used. `origcoords` below is simply a vector of complex numbers containing the original mesh coordinates which will be Delaunay triangulated. For maximum efficiency, the original mesh nodes should form equilateral triangles.
```julia
using GRPF

origcoords = rectangulardomain(complex(xb, yb), complex(xe, ye), r)
```

Roots and poles can be obtained with the `grpf` function. We only need to pass the handle to our `simplefcn` and the `origcoords`.
```julia
zroots, zpoles = grpf(simplefcn, origcoords)
```

### Additional parameters

Additional parameters can be provided to the tesselation and GRPF algorithms by explicitly passing a `GRPFParams` struct.

The two most useful parameters are `tess_sizehint` for the final total number of nodes in the internal `DelaunayTessellation2D` object and the root finder `tolerance` at which the mesh refinement stops. Just like `sizehint!` for other collections, setting `tess_sizehint` to a value approximately equal to the final number of nodes in the tesselation can improve performance. `tolerance` is the largest triangle edge length of the candidate edges defined in the `origcoords` domain. In practice, the root and pole accuracy is a larger value than `tolerance`, so `tolerance` needs to be set smaller than the desired tolerance on the roots and poles.

By default, the value of `tess_sizehint` is 5000 and the `tolerance` is 1e-9, but they can be specified by providing the `GRPFParams` argument
```julia
zroots, zpoles = grpf(simplefcn, origcoords, GRPFParams(5000, 1e-9))
```

Additional parameters which can be controlled are `maxiterations`, `maxnodes`, and `skinnytriangle`. `maxiterations` sets the maximum number of mesh refinement iterations and `maxnodes` sets the maximum number of nodes allowed in the `DelaunayTessellation2D` before returning. `skinnytriangle` is the maximum allowed ratio of the longest to shortest side length in a tesselation triangle before the triangle is automatically subdivided in the mesh refinement step. Default values are

  - `maxiterations`: 100
  - `maxnodes`: 500000
  - `skinnytriangle`: 3

These can be specified along with the `tess_sizehint` and `tolerance` as, e.g.
```julia
zroots, zpoles = grpf(simplefcn, origcoords, tolerance, GRPFParams(100, 500000, 3, 5000, 1e-9))
```

### Plot data

If mesh node `quadrants` and `phasediffs` are wanted for plotting, simply pass a `PlotData()` instance.
```julia
zroots, zpoles, quadrants, phasediffs = grpf(graphenefunction, origcoords, PlotData())
```

A `GRPFParams` can also be passed.
```julia
zroots, zpoles, quadrants, phasediffs = grpf(graphenefunction, origcoords, PlotData(), GRPFParams(5000, 1e-9))
```

### Additional examples

See [test/](test/) for additional examples.

## Limitations

This package uses [VoronoiDelaunay.jl](https://github.com/JuliaGeometry/VoronoiDelaunay.jl) to perform the Delaunay tesselation. `VoronoiDelaunay` is numerically limited to the range of `1.0+eps(Float64)` to `2.0-2eps(Float64)` for its point coordinates. `GRPF.jl` will accept functions and `origcoords` that aren't limited to `Complex{Float64}`, for example `Complex{BigFloat}`, but the internal tolerance of the root finding is limited to `Float64` precision.

## Citing

Please consider citing Piotr's publications if this code is used in scientific work:

  1. P. Kowalczyk, “Complex Root Finding Algorithm Based on Delaunay Triangulation”, ACM Transactions on Mathematical Software, vol. 41, no. 3, art. 19, pp. 1-13, June 2015. https://dl.acm.org/citation.cfm?id=2699457

  2. P. Kowalczyk, "Global Complex Roots and Poles Finding Algorithm Based on Phase Analysis for Propagation and Radiation Problems," IEEE Transactions on Antennas and Propagation, vol. 66, no. 12, pp. 7198-7205, Dec. 2018. https://ieeexplore.ieee.org/document/8457320

We also encourage you to cite this package if used in scientific work. Refer to the Zenodo DOI at the top of the page or [CITATION.bib](CITATION.bib).
