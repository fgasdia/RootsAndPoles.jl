function simplefcn(z)
    w = (z - 1)*(z - im)^2*(z + 1)^3/(z + im)
end

true_zroots = [-1.0 - 0.0im,
                0.0 + 1.0im,
                1.0 - 0.0im]

true_zpoles = [0.0 - 1.0im]

origcoords = ComplexMesh([-2-2im, 2-2im, -2+2im, 2+2im]; rng=RNG)

coords = deepcopy(origcoords)
zroots, zpoles = rootsandpoles(simplefcn, coords)

@test approxmatch(zroots, true_zroots)
@test approxmatch(zpoles, true_zpoles)

# Test a region without any roots or poles
coords = ComplexMesh([-0.2+0.8im, 0.2+0.8im, -0.2+0.9im, 0.2+0.9im]; rng=RNG)
zroots, zpoles = rootsandpoles(simplefcn, coords)

@test length(zroots) == 0
@test length(zpoles) == 0
