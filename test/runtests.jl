using Test
using LinearAlgebra, Random
using SpecialFunctions, StableRNGs
using StaticArrays
using DelaunayTriangulation
const DT = DelaunayTriangulation

using RootsAndPoles
const RP = RootsAndPoles

const RNG = StableRNG(123)

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
    # include("RootsAndPoles.jl")

    @testset "Simple Rational Function" begin include("SimpleRationalFunction.jl") end
    @testset "Default Function" begin include("DefaultFunction.jl") end
    @testset "Complex Modes" begin include("ComplexModes.jl") end
    @testset "Lossy Multilayered Waveguide" begin include("LossyMultilayeredWaveguide.jl") end
    @testset "Graphene Transmission Line" begin include("GrapheneTransmissionLine.jl") end
end
