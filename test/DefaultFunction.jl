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

rmin, rmax = minimum(real(origcoords)), maximum(real(origcoords))
imin, imax = minimum(imag(origcoords)), maximum(imag(origcoords))

ra = (max_coord-min_coord)/(rmax-rmin)
rb = max_coord - ra*rmax

ia = (max_coord-min_coord)/(imax-imin)
ib = max_coord - ia*imax

origcoords = mapfunctionval.(origcoords, ra, rb, ia, ib)
newnodes = [IndexablePoint2D(real(coord), imag(coord), idx) for (idx, coord) in enumerate(origcoords)]
tess = DelaunayTessellation2D{IndexablePoint2D}(2000)

tess, 𝓔, quadrants = GRPF.tesselate!(tess, newnodes, pt -> defaultfcn(geom2fcn(pt, ra, rb, ia, ib)),
                                   e -> geom2fcn(e, ra, rb, ia, ib), tolerance)

𝐶 = GRPF.contouredges(tess, 𝓔)
regions = GRPF.evaluateregions!(𝐶, e -> geom2fcn(e, ra, rb, ia, ib))

zroots, zpoles = GRPF.rootsandpoles(regions, quadrants, e -> geom2fcn(e, ra, rb, ia, ib))

@test length(zroots) == 6
@test length(zpoles) == 2

matlab_zroots = [-1.624715288135189 + 0.182095877702038im,
                 -1.520192978034417 - 0.173670452237129im,
                 -0.515113098919392 + 0.507111597359180im,
                  0.515113098795215 - 0.507111597284675im,
                  1.520192978034417 + 0.173670452237129im,
                  1.624715288135189 - 0.182095877702037im]

matlab_zpoles = [-1.570796326699632 - 0.000000000206961im,
                  1.570796326699632 - 0.000000000206961im]

@test approxmatch(zroots, matlab_zroots)
@test approxmatch(zpoles, matlab_zpoles)

# grpf()
newnodes = [IndexablePoint2D(real(coord), imag(coord), idx) for (idx, coord) in enumerate(origcoords)]
tess = DelaunayTessellation2D{IndexablePoint2D}(2000)

gzroots, gzpoles = grpf(tess, newnodes, pt -> defaultfcn(geom2fcn(pt, ra, rb, ia, ib)),
                        e -> geom2fcn(e, ra, rb, ia, ib), tolerance)

@test approxmatch(zroots, gzroots)
@test approxmatch(zpoles, gzpoles)

newnodes = [IndexablePoint2D(real(coord), imag(coord), idx) for (idx, coord) in enumerate(origcoords)]
tess = DelaunayTessellation2D{IndexablePoint2D}(2000)
gzrootspd, gzpolespd = grpf(tess, newnodes, pt -> defaultfcn(geom2fcn(pt, ra, rb, ia, ib)),
                          e -> geom2fcn(e, ra, rb, ia, ib), tolerance, PhaseDiffs())

@test approxmatch(gzroots, gzrootspd)
@test approxmatch(gzpoles, gzpolespd)
