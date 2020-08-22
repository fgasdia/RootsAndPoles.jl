using BenchmarkTools
using LinearAlgebra

using RootsAndPoles

# Define a parent BenchmarkGroup to contain our suite
const suite = BenchmarkGroup()

# Add some child groups to our benchmark suite.
suite["grpf"] = BenchmarkGroup(["functions"])  # BenchmarkGroup "tags"

simplefcn(z) = (z - 1)*(z - im)^2*(z + 1)^3/(z + im)

function defaultfcn(z)
    f = 1e9
    ϵᵣ = 5 - 2im
    μᵣ = 1 - 2im
    d = 1e-2
    c = 3e8
    ω = 2π*f
    k₀ = ω/c
    cc = ϵᵣ^2*(k₀*d)^2*(ϵᵣ*μᵣ - 1)
    w = ϵᵣ^2*z^2 + z^2*tan(z)^2 - cc
end

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
           im*γs  -m21+im*γc*m22])
end

function graphenefunction(z)
    f = 1e12
    c = 299792458
    μ₀ = 4π*1e-7
    ϵ₀ = 1/(μ₀*c^2)

    e = 1.602176565e-19
    kB = 1.3806488e-23
    hk = 1.05457168e-34
    vFe = 1e6
    muc = 0.05*e
    t = 0.135e-12
    T = 300
    ϵᵣ₁ = 1.0
    ϵᵣ₂ = 11.9

    ω = 2π*f
    k₀ = ω/c
    kᵣ₀ = -im*z*k₀

    Slo=-im*e^2*kB*T*log(2+2*cosh(muc/kB/T)) / (π*hk^2*(ω-im/t))

    a = -3*vFe^2*Slo/(4*(ω-im/t)^2)
    b = a/3

    Y1TM = ω*ϵᵣ₁*ϵ₀/sqrt(ϵᵣ₁*k₀^2 - kᵣ₀^2);
    Y2TM = ω*ϵᵣ₂*ϵ₀/sqrt(ϵᵣ₂*k₀^2 - kᵣ₀^2);
    YSTM = Slo + 1*a*kᵣ₀^2 + 1*b*kᵣ₀^2;

    w = (Y1TM + Y2TM + YSTM)*(-Y1TM + Y2TM + YSTM)*(Y1TM - Y2TM + YSTM)*(-Y1TM - Y2TM + YSTM) # four Riemann sheets
end

testfcns = (simplefcn, defaultfcn, wvgd, graphenefunction)

# Add some benchmarks to the "grpf" group
for f in testfcns
    sf = string(f)
    if sf == "simplefcn"
        origcoords = rectangulardomain(complex(-2, -2), complex(2, 2), 0.1)
        suite["grpf"][sf] = @benchmarkable grpf($f, $origcoords)
    elseif sf == "defaultfcn"
        origcoords = rectangulardomain(complex(-2, -2), complex(2, 2), 0.2)
        suite["grpf"][sf] = @benchmarkable grpf($f, $origcoords)
    elseif sf == "wvgd"
        origcoords = rectangulardomain(complex(1.0, -1.0), complex(2.5, 2.5), 0.5)
        suite["grpf"][sf] = @benchmarkable grpf($f, $origcoords)
    elseif sf == "graphenefunction"
        origcoords = rectangulardomain(complex(-100, -100), complex(400, 400), 18)
        suite["grpf"][sf] = @benchmarkable grpf($f, $origcoords)
    end
end

# If a cache of tuned parameters already exists, use it, otherwise, tune and cache
# the benchmark parameters. Reusing cached parameters is faster and more reliable
# than re-tuning `suite` every time the file is included.
paramspath = joinpath(dirname(@__FILE__), "params.json")

if isfile(paramspath)
    loadparams!(suite, BenchmarkTools.load(paramspath)[1], :evals);
else
    tune!(suite)
    BenchmarkTools.save(paramspath, params(suite));
end

# run with a time limit of ~5 second per benchmark
results = run(suite, verbose = true, seconds = 5)

# Appends results into `paramspath`
BenchmarkTools.save(paramspath, results)
