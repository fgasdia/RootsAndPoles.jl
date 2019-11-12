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

ggzroots, ggzpoles = grpf(defaultfcn, origcoords, tolerance)

@test approxmatch(ggzroots, matlab_zroots)
@test approxmatch(ggzpoles, matlab_zpoles)

ggpzroots, ggpzpoles, quadrants, phasediffs = grpf(defaultfcn, origcoords, tolerance, PlotData())

@test approxmatch(ggpzroots, matlab_zroots)
@test approxmatch(ggpzpoles, matlab_zpoles)


#==
More specific tests
==#
rmin, rmax = minimum(real(origcoords)), maximum(real(origcoords))
imin, imax = minimum(imag(origcoords)), maximum(imag(origcoords))

ra = (max_coord-min_coord)/(rmax-rmin)
rb = max_coord - ra*rmax

ia = (max_coord-min_coord)/(imax-imin)
ib = max_coord - ia*imax

origcoords = GRPF.fcn2geom.(origcoords, ra, rb, ia, ib)
newnodes = [IndexablePoint2D(real(coord), imag(coord), idx) for (idx, coord) in enumerate(origcoords)]
tess = DelaunayTessellation2D{IndexablePoint2D}(2000)

f = GRPF.ScaledFunction(defaultfcn, ra, rb, ia, ib)
g2f = GRPF.Geometry2Function(ra, rb, ia, ib)

tess, 𝓔, quadrants = GRPF.tesselate!(tess, newnodes, f, g2f, tolerance)

𝐶 = GRPF.contouredges(tess, 𝓔)
regions = GRPF.evaluateregions!(𝐶, g2f)

zroots, zpoles = GRPF.rootsandpoles(regions, quadrants, g2f)

@test length(zroots) == 6
@test length(zpoles) == 2

@test approxmatch(zroots, matlab_zroots)
@test approxmatch(zpoles, matlab_zpoles)
