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

# matlab results from https://github.com/PioKow/GRPF for comparison
matlab_zroots = [1.096752543421462 - 0.000197146739811im,
                 1.240454471623525 - 0.000133821833879im,
                 1.353140429314226 - 0.000086139142513im,
                 1.439795544324443 - 0.000052001606673im,
                 1.504169866163284 - 0.000028029270470im,
                 1.548692244058475 - 0.000012100953609im,
                 1.574863045662642 - 0.000002974364907im]

matlab_zpoles = ComplexF64[]

ggzroots, ggzpoles = grpf(wvgd, origcoords, tolerance)

@test approxmatch(ggzroots, matlab_zroots)
@test approxmatch(ggzpoles, matlab_zpoles)

ggpzroots, ggpzpoles, quadrants, phasediffs, tess = grpf(wvgd, origcoords, tolerance, PlotData())

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

f = GRPF.ScaledFunction(wvgd, ra, rb, ia, ib)
g2f = GRPF.Geometry2Function(ra, rb, ia, ib)

tess, 𝓔, quadrants = GRPF.tesselate!(tess, newnodes, f, g2f, tolerance)

𝐶 = GRPF.contouredges(tess, 𝓔)
regions = GRPF.evaluateregions!(𝐶, g2f)

zroots, zpoles = GRPF.rootsandpoles(regions, quadrants, g2f)

@test length(zroots) == 7
@test length(zpoles) == 0

@test approxmatch(zroots, matlab_zroots)
@test approxmatch(zpoles, matlab_zpoles)
