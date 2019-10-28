function wvgd(z)
      ns = 0.065-4im
      n1 = 1.5835
      nc = 1.0
      d1 = 1.81e-6
      λ₀ = 0.6328e-6
      k₀ = 2π/λ₀
      k₀d1 = k₀*d1
      κ1 = sqrt(n1^2 - z^2)
      γs = sqrt(z^2 - ns^2)
      γc = sqrt(z^2 - nc^2)
      m11 = cos(κ1*k₀d1)
      m12 = im/κ1*sin(κ1*k₀d1)
      m21 = im*κ1*sin(κ1*k₀d1)
      m22 = cos(κ1*k₀d1)
      w = det([1.0    -m11+im*γc*m12
               im*γs  -m21+im*γc*m22])
end

# Analysis parameters
xb = 1.0  # real part begin
xe = 2.5  # real part end
yb = -1.0  # imag part begin
ye = 1.0  # imag part end
r = 0.5  # initial mesh step
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

tess, 𝓔, quadrants = GRPF.tesselate!(tess, newnodes, pt -> wvgd(geom2fcn(pt, ra, rb, ia, ib)),
                                     e -> geom2fcn(e, ra, rb, ia, ib), tolerance)

𝐶 = GRPF.contouredges(tess, 𝓔)
regions = GRPF.evaluateregions!(𝐶, e -> geom2fcn(e, ra, rb, ia, ib))

zroots, zpoles = GRPF.rootsandpoles(regions, quadrants, e -> geom2fcn(e, ra, rb, ia, ib))

sort!(zroots, by = x -> (real(x), imag(x)))
sort!(zpoles, by = x -> (real(x), imag(x)))

@test length(zroots) == 7
@test length(zpoles) == 0

@test zroots[1] ≈ 1.096752543421462 - 0.000197146739811im
@test zroots[2] ≈ 1.240454471623525 - 0.000133821833879im
@test zroots[3] ≈ 1.353140429314226 - 0.000086139142513im
@test zroots[4] ≈ 1.439795544324443 - 0.000052001606673im
@test zroots[5] ≈ 1.504169866163284 - 0.000028029270470im
@test zroots[6] ≈ 1.548692244058475 - 0.000012100953609im
@test zroots[7] ≈ 1.574863045662642 - 0.000002974364907im
