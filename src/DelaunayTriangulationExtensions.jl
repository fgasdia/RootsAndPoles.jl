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
setquadrant!(p::QuadrantPoint, q) = (p.q = q)
Base.complex(p::QuadrantPoint) = complex(p.x, p.y)

DT.getx(p::QuadrantPoint) = p.x
DT.gety(p::QuadrantPoint) = p.y
DT.number_type(::Type{QuadrantPoint}) = Float64

DT.each_point_index(pts::QuadrantPoints) = eachindex(pts.points)
DT.num_points(pts::QuadrantPoints) = length(pts.points)
DT.each_point(pts::QuadrantPoints) = pts.points
DT.number_type(::Type{QuadrantPoints}) = DT.number_type(QuadrantPoint)
DT.getpoint(pts::QuadrantPoints, i::Integer) = pts.points[i]
DT.push_point!(pts::QuadrantPoints, pt::QuadrantPoint) = push!(pts.points, pt)
Base.iterate(pts::QuadrantPoints, state...) = Base.iterate(pts.points, state...)
Base.length(pts::QuadrantPoints) = Base.length(pts.points)
Base.firstindex(pts::QuadrantPoints) = 1
Base.lastindex(pts::QuadrantPoints) = length(pts)
Base.getindex(pts::QuadrantPoints, i::Integer) = pts.points[i]
Base.empty!(pts::QuadrantPoints) = empty!(pts.points)
