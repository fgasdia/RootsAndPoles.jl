function complexmodes(z)
      z *= 10
      f = 5e9
      c = 3e8
      μ₀ = 4e-7π
      ϵ₀ = 1e-9/36/π
      a = 6.35e-3
      b = 10e-3
      ϵᵣ₁ = 10
      ϵᵣ₂ = 1
      m = 1

      ω = 2π*f
      k₀ = ω/c
      α = real(z)*k₀
      β = imag(z)*k₀
      γ = α + im*β
      ϵ₁ = ϵ₀*ϵᵣ₁
      ϵ₂ = ϵ₀*ϵᵣ₂
      μ₁ = μ₀
      μ₂ = μ₀
      κ₁ = sqrt(γ^2 + k₀^2*ϵᵣ₁)
      κ₂ = sqrt(γ^2 + k₀^2*ϵᵣ₂)
      η₁ = sqrt(μ₁/ϵ₁)
      η₂ = sqrt(μ₂/ϵ₂)

      Jm_a1 = besselj1(κ₁*a)
      Jm_a2 = besselj1(κ₂*a)
      Ym_a2 = bessely1(κ₂*a)
      Jm_b2 = besselj1(κ₂*b)
      Ym_b2 = bessely1(κ₂*b)
      DJm_a1 = (besselj0(κ₁*a) - besselj(2, κ₁*a))/2
      DJm_a2 = (besselj0(κ₂*a) - besselj(2, κ₂*a))/2
      DJm_b2 = (besselj0(κ₂*b) - besselj(2, κ₂*b))/2
      DYm_a2 = (bessely0(κ₂*a) - bessely(2, κ₂*a))/2
      DYm_b2 = (bessely0(κ₂*b) - bessely(2, κ₂*b))/2

      W = [Jm_a1                 0                       -Jm_a2              -Ym_a2              0                       0;
           0                     Jm_a1/η₁                0                   0                   -Jm_a2/η₂               -Ym_a2/η₂;
           γ*m*Jm_a1/(a*κ₁^2)    -ω*μ₁*DJm_a1/(κ₁*η₁)    -γ*m*Jm_a2/(a*κ₂^2) -γ*m*Ym_a2/(a*κ₂^2) ω*μ₂*DJm_a2/(κ₂*η₂)     ω*μ₂*DYm_a2/(κ₂*η₂);
           -ω*ϵ₁*DJm_a1/κ₁       -m*γ*Jm_a1/(a*κ₁^2*η₁)  ω*ϵ₂*DJm_a2/κ₂      ω*ϵ₂*DYm_a2/κ₂      m*γ*Jm_a2/(a*κ₂^2*η₂)   m*γ*Ym_a2/(a*κ₂^2*η₂);
           0                     0                       Jm_b2               Ym_b2               0                       0;
           0                     0                       γ*m*Jm_b2/(b*κ₂^2)  γ*m*Ym_b2/(b*κ₂^2)  -ω*μ₂*DJm_b2/(κ₂*η₂)    -ω*μ₂*DYm_b2/(κ₂*η₂)]
      w = det(W)
end

R = 1.0
r = 0.15
tolerance = 1e-9

origcoords = diskdomain(R, r)

# matlab results from https://github.com/PioKow/GRPF for comparison
matlab_zroots = [-0.856115203791905 + 0.000000000114004im,
                 -0.775021521974263 + 0.000000000017321im,
                 -0.703772250164549 - 0.000000000209979im,
                 -0.444429043571261 - 0.000000000174298im,
                 -0.096642302493522 - 0.062923397621120im,
                 -0.096642302349294 + 0.062923397246964im,
                  0.096642302349294 - 0.062923397246964im,
                  0.096642302363736 + 0.062923397315544im,
                  0.444429043571261 + 0.000000000174298im,
                  0.703772250003777 - 0.000000000078642im,
                  0.775021521891359 - 0.000000000075087im,
                  0.856115203873180 - 0.000000000219841im]

matlab_zpoles = [0.000000000022863 - 0.100000000307523im,
                 0.000000000069120 + 0.100000000307523im]

ggzroots, ggzpoles = grpf(complexmodes, origcoords, tolerance)

@test approxmatch(ggzroots, matlab_zroots)
@test approxmatch(ggzpoles, matlab_zpoles)

ggpzroots, ggpzpoles, quadrants, phasediffs = grpf(complexmodes, origcoords, tolerance, PlotData())

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

f = GRPF.ScaledFunction(complexmodes, ra, rb, ia, ib)
g2f = GRPF.Geometry2Function(ra, rb, ia, ib)
tess, 𝓔, quadrants = GRPF.tesselate!(tess, newnodes, f, g2f, tolerance)

𝐶 = GRPF.contouredges(tess, 𝓔)
regions = GRPF.evaluateregions!(𝐶, g2f)

zroots, zpoles = GRPF.rootsandpoles(regions, quadrants, g2f)

@test length(zroots) == 12
@test length(zpoles) == 2

@test approxmatch(zroots, matlab_zroots)
@test approxmatch(zpoles, matlab_zpoles)
