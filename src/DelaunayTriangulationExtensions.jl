mutable struct QuadrantPoint{T<:AbstractFloat}
    const x::T
    const y::T
    q::Int
end
QuadrantPoint(x, y) = QuadrantPoint(float(x), float(y), 0)
QuadrantPoint(xy::Tuple) = QuadrantPoint(xy...)
QuadrantPoint(xy::Complex) = QuadrantPoint(real(xy), imag(xy))

struct QuadrantPoints{T}
    points::Vector{QuadrantPoint{T}}
end

getq(p::QuadrantPoint) = p.q
setquadrant!(p::QuadrantPoint, q) = (p.q = q)
get_quadrant(pts::QuadrantPoints, i) = getq(pts.points[i])

get_qpoint(pts::QuadrantPoints, i::Integer) = pts.points[i]  # possibly a bug/improvement, but DT hard codes tuple pairs (x,y) to check things when calling get_point

DT.getx(p::QuadrantPoint) = p.x
DT.gety(p::QuadrantPoint) = p.y
DT.number_type(::Type{QuadrantPoint{T}}) where T = T

Base.eachindex(pts::QuadrantPoints) = Base.eachindex(pts.points)
Base.iterate(pts::QuadrantPoints, state...) = Base.iterate(pts.points, state...)
Base.length(pts::QuadrantPoints) = Base.length(pts.points)
Base.getindex(pts::QuadrantPoints, i) = pts.points[i]
Base.firstindex(pts::QuadrantPoints) = 1
Base.lastindex(pts::QuadrantPoints) = length(pts)
Base.empty!(pts::QuadrantPoints) = empty!(pts.points)

DT.number_type(::Type{QuadrantPoints{T}}) where T = T
DT.is_planar(::QuadrantPoints) = true
