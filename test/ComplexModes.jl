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
        -ω*ϵ₁*DJm_a1/κ₁       -m*γ*Jm_a1/(a*κ₁^2*η₁)   ω*ϵ₂*DJm_a2/κ₂      ω*ϵ₂*DYm_a2/κ₂      m*γ*Jm_a2/(a*κ₂^2*η₂)   m*γ*Ym_a2/(a*κ₂^2*η₂);
         0                     0                       Jm_b2               Ym_b2               0                       0;
         0                     0                       γ*m*Jm_b2/(b*κ₂^2)  γ*m*Ym_b2/(b*κ₂^2)  -ω*μ₂*DJm_b2/(κ₂*η₂)    -ω*μ₂*DYm_b2/(κ₂*η₂)]
    w = det(W)
end

# results from https://github.com/PioKow/GRPF for comparison
# for tol=1e-9, 3867 nodes and 30 iterations of regular GRPF
true_zroots = [-0.09664230246 - 0.06292339746im,
               -0.09664230246 + 0.06292339746im,
                0.09664230246 - 0.06292339746im,
                0.09664230246 + 0.06292339746im,
               -0.44442904311 + 0.0im,
                0.44442904311 - 0.0im,
               -0.70377225022 + 0.0im,
                0.70377225022 - 0.0im,
               -0.77502152220 + 0.0im,
                0.77502152220 - 0.0im,
               -0.85611520391 + 0.0im,
                0.85611520391 - 0.0im]

true_zpoles = [0.0 - 0.1im,
               0.0 + 0.1im]

matlab_zrm = ones(length(true_zroots))
matlab_zpm = [-2, -2]



origcoords = ComplexMesh([-1-1im, 1-1im, 1+1im, -1+1im]; rng=RNG)

coords = deepcopy(origcoords)
zroots, zpoles = rootsandpoles(complexmodes, coords)

@test approxmatch(zroots, true_zroots)
@test approxmatch(zpoles, true_zpoles)

coords = deepcopy(origcoords)
iterations = MeshIterations(coords)
zroots, zpoles = rootsandpoles(complexmodes, coords; iterations, params=FinderParams(maxadaptivenodes=3000))

@test approxmatch(zroots, true_zroots)
@test approxmatch(zpoles, true_zpoles)

R = 1.0
r = 0.15
densecoords = ComplexMesh(RP.diskdomain(R, r); rng=RNG)
zroots, zpoles = rootsandpoles(complexmodes, densecoords; params=FinderParams(maxadaptivenodes=1000))

@test approxmatch(zroots, true_zroots)
@test approxmatch(zpoles, true_zpoles)
