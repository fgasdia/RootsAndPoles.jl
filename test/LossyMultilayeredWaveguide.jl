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
             1im*γs  -m21+im*γc*m22])
end

# matlab results from https://github.com/PioKow/GRPF for comparison
# for tol=1e-9, 2472 nodes and 31 iterations for regular GRPF
true_zroots = [1.09675254341 - 0.00019714688im,
               1.24045447136 - 0.00013382215im,
               1.35314042918 - 0.00008613919im,
               1.43979554425 - 0.00005200167im,
               1.50416986640 - 0.00002802944im,
               1.54869224388 - 0.00001210101im,
               1.57486304575 - 0.00000297462im]

true_zpoles = ComplexF64[]

true_zrm = ones(length(true_zroots))
true_zpm = []

origcoords = ComplexMesh([1-1im, 2.5-1im, 1+1im, 2.5+1im]; rng=RNG)

coords = deepcopy(origcoords)
zroots, zpoles = rootsandpoles(wvgd, coords)

@test approxmatch(zroots, true_zroots)
@test approxmatch(zpoles, true_zpoles)

# ReturnMultiplicity
coords = deepcopy(origcoords)
zroots, zpoles, zrm, zpm = rootsandpoles(wvgd, coords, ReturnMultiplicity())

@test approxmatch(zroots, true_zroots)
@test approxmatch(zpoles, true_zpoles)

@test zrm == true_zrm
@test zpm == true_zpm

# MeshIterations
coords = deepcopy(origcoords)
iterations = MeshIterations(coords)
zroots, zpoles = rootsandpoles(wvgd, coords; iterations)

@test approxmatch(zroots, true_zroots)
@test approxmatch(zpoles, true_zpoles)

# multithreaded
params = FinderParams(numtasks=Threads.nthreads())
coords = deepcopy(origcoords)
zroots, zpoles = rootsandpoles(wvgd, coords; params)

@test approxmatch(zroots, true_zroots)
@test approxmatch(zpoles, true_zpoles)
