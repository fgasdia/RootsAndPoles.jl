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

rmin, rmax = minimum(real(origcoords)), maximum(real(origcoords))
imin, imax = minimum(imag(origcoords)), maximum(imag(origcoords))

ra = (max_coord-min_coord)/(rmax-rmin)
rb = max_coord - ra*rmax

ia = (max_coord-min_coord)/(imax-imin)
ib = max_coord - ia*imax

origcoords = mapfunctionval.(origcoords, ra, rb, ia, ib)
newnodes = [IndexablePoint2D(real(coord), imag(coord), idx) for (idx, coord) in enumerate(origcoords)]
tess = DelaunayTessellation2D{IndexablePoint2D}(2000)

tess, 𝓔, quadrants = GRPF.tesselate!(tess, newnodes, pt -> graphenefunction(geom2fcn(pt, ra, rb, ia, ib)),
                                     e -> geom2fcn(e, ra, rb, ia, ib), tolerance)

𝐶 = GRPF.contouredges(tess, 𝓔)
regions = GRPF.evaluateregions!(𝐶, e -> geom2fcn(e, ra, rb, ia, ib))

zroots, zpoles = GRPF.rootsandpoles(regions, quadrants, e -> geom2fcn(e, ra, rb, ia, ib))

sort!(zroots, by = x -> (real(x), imag(x)))
sort!(zpoles, by = x -> (real(x), imag(x)))

@test length(zroots) == 8
@test length(zpoles) == 2

@test zroots[1] ≈ -38.1777253145628 - 32.5295210454247im
@test zroots[2] ≈ -32.1019622517269 - 27.4308619361753im
@test zroots[3] ≈ 32.1019622517269 + 27.4308619360714im
@test zroots[4] ≈ 38.17772531429 + 32.5295210455806im
@test zroots[5] ≈ 332.744888929695 + 282.243079954389im
@test zroots[6] ≈ 336.220287339074 + 285.191091013829im
@test zroots[7] ≈ 368.439467215558 + 312.522078059503im
@test zroots[8] ≈ 371.007570834263 + 314.700407676927im

# BUG: Sometimes one of zpoles is ~+0 even though Matlab calculates them as both ~-0.
# This causes zpoles[1] and [2] to be flipped and test fails.
if imag(zpoles[1]) < 0
    @test zpoles[1] ≈ -2.30871731988513e-10 - 3.44963766202144im
    @test zpoles[2] ≈ -2.65852297441317e-10 + 3.4496376622893im
else
    @test zpoles[1] ≈ -2.65852297441317e-10 + 3.4496376622893im
    @test zpoles[2] ≈ -2.30871731988513e-10 - 3.44963766202144im
end
