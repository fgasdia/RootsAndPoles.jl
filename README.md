# GRPF.jl: Global complex Roots and Poles Finding in Julia

[![Build Status](https://travis-ci.com/EP-Guy/GRPF.jl.svg?token=U9y2eEri8JFAZrWUCrwX&branch=master)](https://travis-ci.com/EP-Guy/GRPF.jl) [![Build status](https://ci.appveyor.com/api/projects/status/gioglmp08jcivc7h?svg=true)](https://ci.appveyor.com/project/EP-Guy/grpf-jl) [![DOI](https://zenodo.org/badge/154031378.svg)](https://zenodo.org/badge/latestdoi/154031378)

A Julia implementation of [GRPF](https://github.com/PioKow/GRPF) by Dr. Piotr Kowalczyk.

## Description

GRPF attempts to **find all the zeros and poles of a (complex) function in a fixed region**. These types of problems are frequently encountered in electromagnetics, but the algorithm can also be used for similar problems in e.g. optics, acoustics, etc.

GRPF first samples the function on a triangular mesh through Delaunay triangulation. Candidate regions to search for roots and poles are determined and the discretized [Cauchy's argument principle](https://en.wikipedia.org/wiki/Argument_principle) is applied _without needing the derivative of the function or integration over the contour_. To improve the accuracy of the results, a self-adaptive mesh refinement occurs inside the identified candidate regions.

![graphenetransmissionline](graphenetransmissionline.svg)

## Usage

### Installation

```julia
]add https://github.com/fgasdia/GRPF.jl
```

### Example Problem

Consider a simple transmission line consisting of a thin graphene layer on a silicone substrate. See Section III. C. of Kowalczyk, 2018 (below) for details.

First, define the single (complex) argument function for which we seek roots and poles. The normalized propagation coefficient `z` for TM modes at frequency `f` can be found from the equation `graphenefunction(z)`.
```julia
function graphenefunction(z)
    f = 1e12
    c = 299792458
    μ₀ = 4π*1e-7
    ϵ₀ = 1/(μ₀*c^2)

    e = 1.602176565e-19
    kB = 1.3806488e-23
    hk = 1.05457168e-34
    vFe = 1e6
    μc = 0.05*e
    τ = 0.135e-12
    T = 300
    ϵᵣ₁ = 1.0
    ϵᵣ₂ = 11.9

    ω = 2π*f
    k₀ = ω/c
    kᵣ₀ = -im*z*k₀

    Slo=-im*e^2*kB*T*log(2+2*cosh(μc/kB/T)) / (π*hk^2*(ω-im/τ))

    a = -3*vFe^2*Slo/(4*(ω-im/τ)^2)
    b = a/3

    Y1TM = ω*ϵᵣ₁*ϵ₀/sqrt(ϵᵣ₁*k₀^2 - kᵣ₀^2);
    Y2TM = ω*ϵᵣ₂*ϵ₀/sqrt(ϵᵣ₂*k₀^2 - kᵣ₀^2);
    YSTM = Slo + 1*a*kᵣ₀^2 + 1*b*kᵣ₀^2;

    w = (Y1TM + Y2TM + YSTM)*(-Y1TM + Y2TM + YSTM)*(Y1TM - Y2TM + YSTM)*(-Y1TM - Y2TM + YSTM) # four Riemann sheets

    return w
end
```

Next, define parameters for the initial grid.
```julia
xb = -100  # real part begin
xe = 400  # real part end
yb = -100  # imag part begin
ye = 400  # imag part end
r = 18  # initial mesh step
tolerance = 1e-9
```

This package includes functions for rectangular and disk shaped domains, but any shape can be used. `origcoords` below is simply a vector of complex numbers containing the original mesh coordinates which will be Delaunay triangulated. For maximum efficiency, the original mesh nodes should form equilateral triangles.
```julia
using GRPF

origcoords = rectangulardomain(complex(xb, yb), complex(xe, ye), r)
```

Roots and poles can be obtained with the `grpf` function. We only need to pass the handle to our `graphenefunction`, the `origcoords`, and a `tolerance` at which we stop mesh refinement. Specifically, `tolerance` is the smallest triangle edge length of the candidate edges.
```julia
zroots, zpoles = grpf(graphenefunction, origcoords, tolerance)
```

If mesh node `quadrants` and `phasediffs` are wanted for plotting, simply pass a `PlotData()` instance.
```julia
zroots, zpoles, quadrants, phasediffs = grpf(graphenefunction, origcoords, tolerance, PlotData())
```

See [test/](test/) for additional examples.

## Citing

Please consider citing Piotr's publications if this code is used in scientific work:

  1. P. Kowalczyk, “Complex Root Finding Algorithm Based on Delaunay Triangulation”, ACM Transactions on Mathematical Software, vol. 41, no. 3, art. 19, pp. 1-13, June 2015. https://dl.acm.org/citation.cfm?id=2699457
  2. P. Kowalczyk, "Global Complex Roots and Poles Finding Algorithm Based on Phase Analysis for Propagation and Radiation Problems," IEEE Transactions on Antennas and Propagation, vol. 66, no. 12, pp. 7198-7205, Dec. 2018. https://ieeexplore.ieee.org/document/8457320

We also encourage you to cite this package if used in scientific work. Refer to the Zenodo DOI at the top of the page or [CITATION.bib](CITATION.bib).
