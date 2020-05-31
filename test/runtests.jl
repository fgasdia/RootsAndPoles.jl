using Test
using LinearAlgebra
using SpecialFunctions
using VoronoiDelaunay

using RootsAndPoles
import RootsAndPoles: IndexablePoint2D

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
    @time @testset "Simple Rational Function" begin include("SimpleRationalFunction.jl") end
    @time @testset "Complex Modes" begin include("ComplexModes.jl") end
    @time @testset "Lossy Multilayered Waveguide" begin include("LossyMultilayeredWaveguide.jl") end
    @time @testset "Graphene Transmission Line" begin include("GrapheneTransmissionLine.jl") end
    @time @testset "Default" begin include("DefaultFunction.jl") end
end
