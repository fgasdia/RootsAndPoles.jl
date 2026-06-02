function graphenefunction(z)
    f = 1e12
    c = 299792458.0
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

    # four Riemann sheets
    w = (Y1TM + Y2TM + YSTM)*(-Y1TM + Y2TM + YSTM)*(Y1TM - Y2TM + YSTM)*(-Y1TM - Y2TM + YSTM)
end

# results from SAGRPF paper, single order
true_zroots = [-32.101962251 - 27.430861936im,
                32.101962251 + 27.430861936im,
               -38.177725314 - 32.529521045im,
                38.177725314 + 32.529521045im,  # typo in SAGRPF paper on sign of Im
               332.744888929 + 282.243079954im,
               336.220287339 + 285.191091014im,
               368.439467216 + 312.522078059im,
               371.007570834 + 314.700407677im,
                −0.003206780 − 0.964810358im,
                 0.003206780 + 0.964810358im,
                −0.004526720 + 0.955901829im,
                 0.004526720 - 0.955901829im]

# second order
true_zpoles = [0.000000000 - 1.000000000im,
               0.000000000 + 1.000000000im,
               0.000000000 - 3.449637662im,
               0.000000000 + 3.449637662im]

true_zrm = ones(length(true_zroots))
true_zpm = [-2, -2, -2, -2]

origcoords = ComplexMesh([-100-100im, 400-100im, -100+400im, 400+400im]; rng=RNG)
params = FinderParams(maxadaptivenodes=10000)

coords = deepcopy(origcoords)
zroots, zpoles = rootsandpoles(graphenefunction, coords; params)

@test approxmatch(zroots, true_zroots)
@test approxmatch(zpoles, true_zpoles)

# MeshIterations
coords = deepcopy(origcoords)
iterations = MeshIterations(coords)
zroots, zpoles = rootsandpoles(graphenefunction, coords; iterations, params)

@test approxmatch(zroots, true_zroots)
@test approxmatch(zpoles, true_zpoles)
