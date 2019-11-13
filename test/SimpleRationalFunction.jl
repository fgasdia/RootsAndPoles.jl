function simplefcn(z)
      w = (z - 1)*(z - im)^2*(z + 1)^3/(z + im)
end

# Analysis parameters
xb = -2  # real part begin
xe = 2  # real part end
yb = -2  # imag part begin
ye = 2  # imag part end
r = 0.1  # initial mesh step
tolerance = 1e-9

origcoords = rectangulardomain(complex(xb, yb), complex(xe, ye), r)

# matlab results from https://github.com/PioKow/GRPF for comparison
matlab_zroots = [-0.999999999951224 - 0.000000000028656im,
                  0.000000000253637 + 1.000000000074506im,
                  1.000000000317046 - 0.000000000062088im]

matlab_zpoles = [0.000000000380455 - 0.999999999701977im]

ggzroots, ggzpoles = grpf(simplefcn, origcoords, tolerance)

@test approxmatch(ggzroots, matlab_zroots)
@test approxmatch(ggzpoles, matlab_zpoles)

ggpzroots, ggpzpoles, quadrants, phasediffs, tess = grpf(simplefcn, origcoords, tolerance, PlotData())

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
tess = DelaunayTessellation2D{IndexablePoint2D}(5000)

f = GRPF.ScaledFunction(simplefcn, ra, rb, ia, ib)
g2f = GRPF.Geometry2Function(ra, rb, ia, ib)

tess, 𝓔, quadrants = GRPF.tesselate!(tess, newnodes, f, g2f, tolerance)

𝐶 = GRPF.contouredges(tess, 𝓔)
regions = GRPF.evaluateregions!(𝐶, g2f)

zroots, zpoles = GRPF.rootsandpoles(regions, quadrants, g2f)

@test length(zroots) == 3
@test length(zpoles) == 1

@test approxmatch(zroots, matlab_zroots)
@test approxmatch(zpoles, matlab_zpoles)
