using Test
using LinearAlgebra
using SpecialFunctions
using DelaunayTriangulation

using RootsAndPoles
const RP = RootsAndPoles

function approxmatch(A, B)
    length(A) == length(B) || return false

    matches = Dict{eltype(A),Vector{eltype(B)}}()

    for a in A
        for b in B
            if b ≈ a
                if haskey(matches, a)
                    push!(matches[a], b)
                else
                    matches[a] = [b]
                end
            end
        end
    end

    length(matches) == length(A) || return false
    any(x->length(x)>1, values(matches)) && return false

    return true
end

@testset "RootsAndPoles" begin
    include("RootsAndPoles.jl")

    # @testset "Simple Rational Function" begin include("SimpleRationalFunction.jl") end
    @testset "Complex Modes" begin include("ComplexModes.jl") end
    @testset "Lossy Multilayered Waveguide" begin include("LossyMultilayeredWaveguide.jl") end
    @testset "Graphene Transmission Line" begin include("GrapheneTransmissionLine.jl") end
    @testset "Default" begin include("DefaultFunction.jl") end
end
