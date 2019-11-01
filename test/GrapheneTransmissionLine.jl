function graphenefunction(z)
      f = 1e12
      c = 299792458
      μ₀ = 4π*1e-7
      ϵ₀ = 1/(μ₀*c^2)

      e = 1.602176565e-19
      kB = 1.3806488e-23
      hk = 1.05457168e-34
      vFe = 1e6
      muc = 0.05*e
      t = 0.135e-12
      T = 300
      ϵᵣ₁ = 1.0
      ϵᵣ₂ = 11.9

      ω = 2π*f
      k₀ = ω/c
      kᵣ₀ = -im*z*k₀

      Slo=-im*e^2*kB*T*log(2+2*cosh(muc/kB/T)) / (π*hk^2*(ω-im/t))

      a = -3*vFe^2*Slo/(4*(ω-im/t)^2)
      b = a/3

      Y1TM = ω*ϵᵣ₁*ϵ₀/sqrt(ϵᵣ₁*k₀^2 - kᵣ₀^2);
      Y2TM = ω*ϵᵣ₂*ϵ₀/sqrt(ϵᵣ₂*k₀^2 - kᵣ₀^2);
      YSTM = Slo + 1*a*kᵣ₀^2 + 1*b*kᵣ₀^2;

      w = (Y1TM + Y2TM + YSTM)*(-Y1TM + Y2TM + YSTM)*(Y1TM - Y2TM + YSTM)*(-Y1TM - Y2TM + YSTM) # four Riemann sheets
end

# Analysis parameters
xb = -100  # real part begin
xe = 400  # real part end
yb = -100  # imag part begin
ye = 400  # imag part end
r = 18  # initial mesh step
tolerance = 1e-9

origcoords = rectangulardomain(complex(xb, yb), complex(xe, ye), r)

# matlab results from https://github.com/PioKow/GRPF for comparison
matlab_zroots = [-38.1777253145628 - 32.5295210454247im,
                 -32.1019622517269 - 27.4308619361753im,
                  32.1019622517269 + 27.4308619360714im,
                  38.17772531429 + 32.5295210455806im,
                  332.744888929695 + 282.243079954389im,
                  336.220287339074 + 285.191091013829im,
                  368.439467215558 + 312.522078059503im,
                  371.007570834263 + 314.700407676927im]

matlab_zpoles = [-2.30871731988513e-10 - 3.44963766202144im,
                 -2.65852297441317e-10 + 3.4496376622893im]

ggzroots, ggzpoles = grpf(graphenefunction, origcoords, tolerance)

@test approxmatch(ggzroots, matlab_zroots)
@test approxmatch(ggzpoles, matlab_zpoles)

ggpzroots, ggpzpoles, quadrants, phasediffs = grpf(graphenefunction, origcoords, tolerance, PlotData())

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

tess, 𝓔, quadrants = GRPF.tesselate!(tess, newnodes, pt -> graphenefunction(GRPF.geom2fcn(pt, ra, rb, ia, ib)),
                                     e -> GRPF.geom2fcn(e, ra, rb, ia, ib), tolerance)

𝐶 = GRPF.contouredges(tess, 𝓔)
regions = GRPF.evaluateregions!(𝐶, e -> GRPF.geom2fcn(e, ra, rb, ia, ib))

zroots, zpoles = GRPF.rootsandpoles(regions, quadrants, e -> GRPF.geom2fcn(e, ra, rb, ia, ib))

@test length(zroots) == 8
@test length(zpoles) == 2

@test approxmatch(zroots, matlab_zroots)
@test approxmatch(zpoles, matlab_zpoles)
