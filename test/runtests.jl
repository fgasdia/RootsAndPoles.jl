using Test
using LinearAlgebra
using SpecialFunctions
using DelaunayTriangulation

using RootsAndPoles
const RP = RootsAndPoles

function approxmatch(A, B)
    length(A) == length(B) || return false

    # Need to make sure we don't match to the same root twice. It's happened...
    Acopy = copy(A)
    Bcopy = copy(B)

    # This could be less than O(n²) but it's not worth the bookkeeping effort
    i = 1
    while length(Acopy) > 0
        amatch = false
        for j in eachindex(Bcopy)
            if Acopy[i] ≈ Bcopy[j]
                amatch = true
                deleteat!(Acopy, i)
                deleteat!(Bcopy, j)
                break
            end
        end
        if !amatch
            return false
        end
        i += 1
    end
    return true
end

@testset "RootsAndPoles" begin
    include("RootsAndPoles.jl")

    @testset "Simple Rational Function" begin include("SimpleRationalFunction.jl") end
    @testset "Complex Modes" begin include("ComplexModes.jl") end
    @testset "Lossy Multilayered Waveguide" begin include("LossyMultilayeredWaveguide.jl") end
    @testset "Graphene Transmission Line" begin include("GrapheneTransmissionLine.jl") end
    @testset "Default" begin include("DefaultFunction.jl") end
end
