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

# `true_zroots` derived from Matlab (https://github.com/PioKow/GRPF) for comparison
true_zroots = [-1.624715288135189 + 0.182095877702038im,
               -1.520192978034417 - 0.173670452237129im,
               -0.515113098919392 + 0.507111597359180im,
                0.515113098795215 - 0.507111597284675im,
                1.520192978034417 + 0.173670452237129im,
                1.624715288135189 - 0.182095877702037im]

true_zpoles = [-π/2 - 0.0im,
                π/2 - 0.0im]

origcoords = ComplexMesh([-2-2im, 2-2im, 2+2im, -2+2im]; rng=RNG)

coords = deepcopy(origcoords)
zroots, zpoles = rootsandpoles(defaultfcn, coords)

@test approxmatch(zroots, true_zroots)
@test approxmatch(zpoles, true_zpoles)

# ReturnMultiplicity
coords = deepcopy(origcoords)
zroots, zpoles, zrm, zpm = rootsandpoles(defaultfcn, coords, ReturnMultiplicity())

@test approxmatch(zroots, true_zroots)
@test approxmatch(zpoles, true_zpoles)

# MeshIterations
coords = deepcopy(origcoords)
iterations = MeshIterations(coords)
zroots, zpoles = rootsandpoles(defaultfcn, coords; iterations)

@test approxmatch(zroots, true_zroots)
@test approxmatch(zpoles, true_zpoles)
