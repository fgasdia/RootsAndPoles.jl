using DelaunayTriangulation
using BenchmarkTools
# using Makie

function rectangulardomain(rzbl, izbl, rztr, iztr, Δr)
    X = rztr - rzbl
    Y = iztr - izbl

    n = ceil(Int, Y/Δr)
    dy = Y/n

    ## dx = sqrt(Δr² - (dy/2)²), solved for equilateral triangle
    m = ceil(Int, X/sqrt(Δr^2 - dy^2/4))
    dx = X/m
    half_dx = dx/2

    v = Vector{Tuple{Float64, Float64}}()
    sizehint!(v, (m+1)*(n+1))

    shift = false  # we will displace every other line by dx/2
    for j = 0:n
        y = izbl + dy*j

        for i = 0:m
            x = rzbl + dx*i

            if shift && i == 0
                continue  # otherwise, we shift out of left bound
            elseif shift
                x -= half_dx
            end

            if i == m
                shift = !shift
            end

            push!(v, (x, y))
        end
    end

    return v
end

const DT = DelaunayTriangulation

mutable struct QuadrantPoint
    const x::Float64
    const y::Float64
    q::Int
end
QuadrantPoint(x, y) = QuadrantPoint(x, y, 0)
QuadrantPoint(xy::Tuple) = QuadrantPoint(xy...)
QuadrantPoint(xy::Complex) = QuadrantPoint(real(xy), imag(xy), 0)

struct QuadrantPoints
    points::Vector{QuadrantPoint}
end

getquadrant(p::QuadrantPoint) = p.q

DT.getx(p::QuadrantPoint) = p.x
DT.gety(p::QuadrantPoint) = p.y
DT.number_type(::Type{QuadrantPoint}) = Float64

DT.each_point_index(pts::QuadrantPoints) = eachindex(pts.points)
DT.num_points(pts::QuadrantPoints) = length(pts.points)
DT.each_point(pts::QuadrantPoints) = pts.points
DT.number_type(::Type{QuadrantPoints}) = DT.number_type(QuadrantPoint)
DT.getpoint(pts::QuadrantPoints, i::Integer) = pts.points[i]

# Generate a large triangulation
coords = rectangulardomain(0.1, 0.1, 0.9, 0.9, 0.01)
quadcoords = QuadrantPoints(QuadrantPoint.(coords))

dt = triangulate(quadcoords)
get_points(dt)

add_point!(dt, (0.23, 0.24))
get_points(dt)

cnt = 0
for T in each_triangle(dt)
    i, j, k = indices(T)
    p, q, r = get_point(dt, i, j, k)
    println(p, q, r)
    break
end 

bt, cdata = basictriangulation(coords) # initialize data structures 
bt, cdata = update!(coords, bt, cdata) # after the points have been changed, may incur allocations
