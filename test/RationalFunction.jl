function rationalfcn(z, ϵ)
    za = 1/2 - sqrt(3)/6*1im
    zb = sqrt(3)/3*1im
    zc = -1/2 - sqrt(3)/6*1im

    w = (z - za)*(z - zb - ϵ)/(z - zc)/(z - zb + ϵ)

    return w
end

za = 1/2 - sqrt(3)/6*1im
zb = sqrt(3)/3*1im
zc = -1/2 - sqrt(3)/6*1im
ϵ = 1e-2

true_roots = [za, zb+ϵ]
true_poles = [zc, zb-ϵ]

mesh = ComplexMesh([-1-1im, 1-1im, 1+1im, -1+1im]; rng=RNG)
roots, poles = @inferred rootsandpoles(z->rationalfcn(z, ϵ), mesh)

@test approxmatch(roots, true_roots)
@test approxmatch(poles, true_poles)

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
