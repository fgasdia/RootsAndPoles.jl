function simplefcn(z)
      w = (z - 1)*(z - im)^2*(z + 1)^3/(z + im)
end

# Analysis parameters
xb = -2.  # real part begin
xe = 2.  # real part end
yb = -2.  # imag part begin
ye = 2.  # imag part end
r = 0.1  # initial mesh step
tolerance = 1e-9

origcoords = rectangulardomain(complex(xb, yb), complex(xe, ye), r)

rmin, rmax = minimum(real(origcoords)), maximum(real(origcoords))
imin, imax = minimum(imag(origcoords)), maximum(imag(origcoords))

ra = (max_coord-min_coord)/(rmax-rmin)
rb = max_coord - ra*rmax

ia = (max_coord-min_coord)/(imax-imin)
ib = max_coord - ia*imax

origcoords = mapfunctionval.(origcoords, ra, rb, ia, ib)
newnodes = [IndexablePoint2D(real(coord), imag(coord), idx) for (idx, coord) in enumerate(origcoords)]
tess = DelaunayTessellation2D{IndexablePoint2D}(5000)

tess, 𝓔, quadrants = GRPF.tesselate!(tess, newnodes, pt -> simplefcn(geom2fcn(pt, ra, rb, ia, ib)),
                                     e -> geom2fcn(e, ra, rb, ia, ib), tolerance)

𝐶 = GRPF.contouredges(tess, 𝓔)
regions = GRPF.evaluateregions!(𝐶, e -> geom2fcn(e, ra, rb, ia, ib))

zroots, zroots_multiplicity, zpoles, zpoles_multiplicity = GRPF.rootsandpoles(regions, quadrants, e -> geom2fcn(e, ra, rb, ia, ib))

sort!(zroots, by = x -> (real(x), imag(x)))

@test length(zroots) == 3
@test length(zpoles) == 1

@test zroots[1] ≈ -0.999999999951224 - 0.000000000028656im
@test zroots[2] ≈ 0.000000000253637 + 1.000000000074506im
@test zroots[3] ≈ 1.000000000317046 - 0.000000000062088im

@test zpoles[1] ≈ 0.000000000380455 - 0.999999999701977im
