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
tolerance = 1e-9

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

zroots, zpoles = grpf(defaultfcn, origcoords, tolerance)

@test length(zroots) == 6
@test length(zpoles) == 2

@test approxmatch(zroots, matlab_zroots)
@test approxmatch(zpoles, matlab_zpoles)

pzroots, pzpoles, quadrants, phasediffs, tess = grpf(defaultfcn, origcoords, tolerance, PlotData())

@test approxmatch(pzroots, matlab_zroots)
@test approxmatch(pzpoles, matlab_zpoles)

# Test with big origcoords
xb = big"-2"  # real part begin
xe = big"2"  # real part end
yb = big"-2"  # imag part begin
ye = big"2"  # imag part end
r = big"0.2"  # initial mesh step
tolerance = 1e-9

origcoords = rectangulardomain(complex(xb, yb), complex(xe, ye), r)

bzroots, bzpoles = grpf(defaultfcn, origcoords, tolerance)

@test all(isa.(bzroots, Complex{BigFloat}))
@test all(isa.(bzpoles, Complex{BigFloat}))

@test approxmatch(bzroots, zroots)
@test approxmatch(bzpoles, zpoles)
