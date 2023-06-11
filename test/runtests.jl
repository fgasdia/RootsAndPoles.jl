using Test
using LinearAlgebra
using SpecialFunctions
using DelaunayTriangulation

using RootsAndPoles
const RP = RootsAndPoles

function approxmatch(A::AbstractArray, B::AbstractArray)
    length(A) == length(B) || return false

    # This could be less than O(n²) but it's not worth the bookkeeping effort
    @inbounds for i in eachindex(A)
        amatch = false
        @inbounds for j in eachindex(B)
            if A[i] ≈ B[j]
                amatch = true
                break
            end
        end
        if !amatch
            return false
        end
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
