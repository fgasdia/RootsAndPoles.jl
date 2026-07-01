function simplefcn(z)
    w = (z - 1)*(z - im)^2*(z + 1)^3/(z + im)
end

# matlab results from https://github.com/PioKow/GRPF for comparison
true_zroots = [-0.999999999951224 - 0.000000000028656im,
                0.000000000253637 + 1.000000000074506im,
                1.000000000317046 - 0.000000000062088im]

true_zpoles = [0.000000000380455 - 0.999999999701977im]

coords = ComplexMesh([-2-2im, 2-2im, -2+2im, 2+2im]; rng=RNG)
zroots, zpoles = rootsandpoles(simplefcn, coords)

@test approxmatch(zroots, true_zroots)
@test approxmatch(zpoles, true_zpoles)
