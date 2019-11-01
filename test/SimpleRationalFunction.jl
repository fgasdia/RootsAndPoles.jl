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

rmin, rmax = minimum(real(origcoords)), maximum(real(origcoords))
imin, imax = minimum(imag(origcoords)), maximum(imag(origcoords))

ra = (max_coord-min_coord)/(rmax-rmin)
rb = max_coord - ra*rmax

ia = (max_coord-min_coord)/(imax-imin)
ib = max_coord - ia*imax

origcoords = GRPF.fcn2geom.(origcoords, ra, rb, ia, ib)
newnodes = [IndexablePoint2D(real(coord), imag(coord), idx) for (idx, coord) in enumerate(origcoords)]
tess = DelaunayTessellation2D{IndexablePoint2D}(5000)

tess, 𝓔, quadrants = GRPF.tesselate!(tess, newnodes, pt -> simplefcn(GRPF.geom2fcn(pt, ra, rb, ia, ib)),
                                     e -> GRPF.geom2fcn(e, ra, rb, ia, ib), tolerance)

𝐶 = GRPF.contouredges(tess, 𝓔)
regions = GRPF.evaluateregions!(𝐶, e -> GRPF.geom2fcn(e, ra, rb, ia, ib))

zroots, zpoles = GRPF.rootsandpoles(regions, quadrants, e -> GRPF.geom2fcn(e, ra, rb, ia, ib))

@test length(zroots) == 3
@test length(zpoles) == 1

matlab_zroots = [-0.999999999951224 - 0.000000000028656im,
                  0.000000000253637 + 1.000000000074506im,
                  1.000000000317046 - 0.000000000062088im]

matlab_zpoles = [0.000000000380455 - 0.999999999701977im]

# grpf()
origcoords = rectangulardomain(complex(xb, yb), complex(xe, ye), r)
ggzroots, ggzpoles = grpf(simplefcn, origcoords, tolerance)

@test approxmatch(ggzroots, matlab_zroots)
@test approxmatch(ggzpoles, matlab_zpoles)

ggpzroots, ggpzpoles, quadrants, phasediffs = grpf(simplefcn, origcoords, tolerance, PlotData())

@test approxmatch(ggpzroots, matlab_zroots)
@test approxmatch(ggpzpoles, matlab_zpoles)
