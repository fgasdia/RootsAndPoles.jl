using Test
using LinearAlgebra
using SpecialFunctions
using VoronoiDelaunay

using GRPF

# TODO: Order output?

"""
    mapfunctionval(z, ra, rb, ia, ib)

Linearly map function values within domain from `min_coord` to `max_coord`.
"""
function mapfunctionval(z, ra, rb, ia, ib)
    zr = ra*real(z) + rb
    zi = ia*imag(z) + ib
    return complex(zr, zi)
end

"""
    geom2fcn(pt, ra, rb, ia, ib)

# TODO: Include this function in package (not as test?)

Linearly map geometry values ∈ {`min_coord`, `max_coord`} to domain bounds.

Note: There are floating point errors when converting back and forth.
"""
function geom2fcn(pt::IndexablePoint2D, ra, rb, ia, ib)
    return complex((getx(pt) - rb)/ra, (gety(pt) - ib)/ia)
end
function geom2fcn(edge::VoronoiDelaunay.DelaunayEdge{IndexablePoint2D}, ra, rb, ia, ib)
    return geom2fcn(geta(edge), ra, rb, ia, ib), geom2fcn(getb(edge), ra, rb, ia, ib)
end

@time @testset "Simple Rational Function" begin include("SimpleRationalFunction.jl") end
@time @testset "Complex Modes" begin include("ComplexModes.jl") end
@time @testset "Lossy Multilayered Waveguide" begin include("LossyMultilayeredWaveguide.jl") end
@time @testset "Graphene Transmission Line" begin include("GrapheneTransmissionLine.jl") end
@time @testset "Default" begin include("DefaultFunction.jl") end
