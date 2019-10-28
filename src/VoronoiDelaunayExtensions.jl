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
getx(p::IndexablePoint2D) = p._x
gety(p::IndexablePoint2D) = p._y
Base.getindex(p::IndexablePoint2D) = p._index
Base.setindex!(p::IndexablePoint2D, v::Int) = setfield!(p, :_index, v)
Base.:+(p1::IndexablePoint2D, p2::IndexablePoint2D) = IndexablePoint2D(getx(p1)+getx(p2), gety(p1)+gety(p2), -1)
Base.:/(p1::IndexablePoint2D, n::Real) = IndexablePoint2D(getx(p1)/n, gety(p1)/n, -1)


# Improved performance of `delaunayedges()`
# see: https://github.com/JuliaGeometry/VoronoiDelaunay.jl/issues/47
function delaunayedges_fast(t::DelaunayTessellation2D)
    result = DelaunayEdge[]
    for ix in 2:t._last_trig_index
        tr = t._trigs[ix]
        isexternal(tr) && continue

        ix_na = tr._neighbour_a
        if (ix_na > ix) || isexternal(t._trigs[ix_na])
            push!(result, DelaunayEdge(getb(tr), getc(tr)))
        end
        ix_nb = tr._neighbour_b
        if (ix_nb > ix) || isexternal(t._trigs[ix_nb])
            push!(result, DelaunayEdge(geta(tr), getc(tr)))
        end
        ix_nc = tr._neighbour_c
        if (ix_nc > ix) || isexternal(t._trigs[ix_nc])
            push!(result, DelaunayEdge(geta(tr), getb(tr)))
        end
    end
    return result
end
