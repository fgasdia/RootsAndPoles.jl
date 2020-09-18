function defaultfcn(z)
    f = 1e9
    ϵᵣ = 5 - 2im
    μᵣ = 1 - 2im
    d = 1e-2
    c = 3e8
    ω = 2π*f
    k₀ = ω/c
    cc = ϵᵣ^2*(k₀*d)^2*(ϵᵣ*μᵣ - 1)
    w = ϵᵣ^2*z^2 + z^2*tan(z)^2 - cc
end

# Analysis parameters
xb = -2  # real part begin
xe = 2  # real part end
yb = -2  # imag part begin
ye = 2  # imag part end
r = 0.2  # initial mesh step

origcoords = rectangulardomain(complex(xb, yb), complex(xe, ye), r)

# matlab results from https://github.com/PioKow/GRPF for comparison
matlab_zroots = [-1.624715288135189 + 0.182095877702038im,
                 -1.520192978034417 - 0.173670452237129im,
                 -0.515113098919392 + 0.507111597359180im,
                  0.515113098795215 - 0.507111597284675im,
                  1.520192978034417 + 0.173670452237129im,
                  1.624715288135189 - 0.182095877702037im]

matlab_zpoles = [-1.570796326699632 - 0.000000000206961im,
                  1.570796326699632 - 0.000000000206961im]

zroots, zpoles = grpf(defaultfcn, origcoords)

@test length(zroots) == 6
@test length(zpoles) == 2

@test approxmatch(zroots, matlab_zroots)
@test approxmatch(zpoles, matlab_zpoles)

pzroots, pzpoles, quadrants, phasediffs, tess, g2f = grpf(defaultfcn, origcoords, PlotData())

@test approxmatch(pzroots, matlab_zroots)
@test approxmatch(pzpoles, matlab_zpoles)

tmpr, tmpp = grpf(defaultfcn, origcoords, GRPFParams(8000, 1e-9))

@test approxmatch(tmpr, zroots)
@test approxmatch(tmpp, zpoles)

tmpr, tmpp, quadrants, phasediffs, tess, g2f = grpf(defaultfcn, origcoords, PlotData(), GRPFParams(8000, 1e-9))

@test approxmatch(tmpr, zroots)
@test approxmatch(tmpp, zpoles)

tmpr, tmpp = grpf(defaultfcn, origcoords, GRPFParams(8000, 1e-9, true))

@test approxmatch(tmpr, zroots)
@test approxmatch(tmpp, zpoles)

# Test with very low maxnodes
maxnodes = 10
@test_logs (:warn,"GRPFParams `tess_sizehint` is greater than `maxnodes`") GRPFParams(100, maxnodes, 3, 8000, 1e-9, false)
params = GRPFParams(100, maxnodes, 3, 8000, 1e-9, false)
# `let` block needed for phasediffs to be defined when we don't exit tesselation early
@test_logs (:warn,"GRPFParams.maxnodes reached") grpf(defaultfcn, origcoords, PlotData(), params)

# Test with big origcoords
xb = big"-2"  # real part begin
xe = big"2"  # real part end
yb = big"-2"  # imag part begin
ye = big"2"  # imag part end
r = big"0.2"  # initial mesh step

origcoords = rectangulardomain(complex(xb, yb), complex(xe, ye), r)

bzroots, bzpoles = grpf(defaultfcn, origcoords)

@test all(isa.(bzroots, Complex{BigFloat}))
@test all(isa.(bzpoles, Complex{BigFloat}))

@test approxmatch(bzroots, zroots)
@test approxmatch(bzpoles, zpoles)

# Test region with no roots or poles
xb, xe = -1, 1
yb, ye = -2, -1
origcoords = rectangulardomain(complex(xb, yb), complex(xe, ye), r)

nozroots, nozpoles = grpf(defaultfcn, origcoords)

@test length(nozroots) == 0
@test length(nozpoles) == 0

nozroots, nozpoles, quadrants, phasediffs, tess, g2f = grpf(defaultfcn, origcoords, PlotData())

@test length(nozroots) == 0
@test length(nozpoles) == 0

xb, xe = big"-1", big"1"
yb, ye = big"-2", big"-1"
r = big"0.2"

origcoords = rectangulardomain(complex(xb, yb), complex(xe, ye), r)

bnozroots, bnozpoles = grpf(defaultfcn, origcoords)

@test eltype(bnozroots) == Complex{BigFloat}
@test eltype(bnozpoles) == Complex{BigFloat}

@test length(bnozroots) == 0
@test length(bnozpoles) == 0
