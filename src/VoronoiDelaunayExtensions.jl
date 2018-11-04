import Base
using VoronoiDelaunay
import VoronoiDelaunay: getx, gety, DelaunayEdge, DelaunayTriangle

export IndexablePoint2D

mutable struct IndexablePoint2D <: AbstractPoint2D
    _x::Float64
    _y::Float64
    _index::Int
end
IndexablePoint2D(x::Float64, y::Float64) = IndexablePoint2D(x, y, -1)
getx(p::IndexablePoint2D)::Float64 = p._x
gety(p::IndexablePoint2D)::Float64 = p._y
Base.getindex(p::IndexablePoint2D)::Int = p._index
Base.setindex!(p::IndexablePoint2D, v::Int) = setfield!(p, :_index, v)
Base.:+(p1::IndexablePoint2D, p2::IndexablePoint2D) = IndexablePoint2D(getx(p1)+getx(p2), gety(p1)+gety(p2), -1)
Base.:/(p1::IndexablePoint2D, n::Real) = IndexablePoint2D(getx(p1)/n, gety(p1)/n, -1)
