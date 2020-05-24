"""
    IndexablePoint2D <: AbstractPoint2D

A special 2D point type, compatible with `VoronoiDelaunay`, that carries an
identifying `_index` field.

By default, `_index` is -1.

See also: `setindex!`, `getindex`
"""
mutable struct IndexablePoint2D <: AbstractPoint2D
    _x::Float64
    _y::Float64
    _index::Int
end
IndexablePoint2D(x::Float64, y::Float64) = IndexablePoint2D(x, y, -1)
getx(p::IndexablePoint2D) = p._x
gety(p::IndexablePoint2D) = p._y
Base.:+(p1::IndexablePoint2D, p2::IndexablePoint2D) = IndexablePoint2D(getx(p1)+getx(p2), gety(p1)+gety(p2), -1)
Base.:/(p1::IndexablePoint2D, n::Real) = IndexablePoint2D(getx(p1)/n, gety(p1)/n, -1)

"""
    getindex(p::IndexablePoint2D)

Return `_index` field of `p`.
"""
Base.getindex(p::IndexablePoint2D) = p._index

"""
    setindex!(p::IndexablePoint2D)

Set the `_index` field of `p`.

# Example

```jldoctest
julia> p = IndexablePoint2D(103.4, 25.0);

julia> getindex(p)
-1

julia> setindex!(p, 90);

julia> getindex(p)
90
```
"""
Base.setindex!(p::IndexablePoint2D, v::Int) = setfield!(p, :_index, v)

"""
    same(e1::DelaunayEdge{IndexablePoint2D}, e2::DelaunayEdge{IndexablePoint2D})

Return `true` if edge `e1` is the same as `e2`, regardless of direction.
Otherwise, return `false`. This is similar to `==`.

To check if two edges are identical in the sense that no program could
distinguish them, use `===`.
"""
function same(e1::DelaunayEdge{IndexablePoint2D},e2::DelaunayEdge{IndexablePoint2D})
    e1 === e2 && return true

    # if they're reversed
    e1a, e1b = geta(e1)._index, getb(e1)._index
    e2a, e2b = geta(e2)._index, getb(e2)._index
    (e1a == e2b) && (e1b == e2a) && return true

    return false
end

"""
    sameunique!

In-place `unique!` that considers both directions of edges the same (just like
[`same`](@ref)).
"""
function sameunique!(v::AbstractVector{DelaunayEdge{IndexablePoint2D}})
    deleteidxs = trues(length(v))
    for i in eachindex(v)
        @inbounds edgei = C[i]
        eia = getindex(geta(edgei))
        eib = getindex(getb(edgei))
        revdupe = false
        for j in eachindex(v)
            @inbounds edgej = v[j]
            eja = getindex(geta(edgej))
            ejb = getindex(getb(edgej))
            if (eia == ejb) && (eib == eja)
                revdupe = true
                break
            end
        end
        if !revdupe
            deleteidxs[i] = false
        end
    end

    deleteat!(v, deleteidxs)
end

# Improved performance of `delaunayedges()`
# see: https://github.com/JuliaGeometry/VoronoiDelaunay.jl/issues/47
function delaunayedges_fast(t::DelaunayTessellation2D{T}) where T <: AbstractPoint2D
    result = DelaunayEdge{T}[]
    @inbounds for ix in 2:t._last_trig_index
        tr = t._trigs[ix]
        isexternal(tr) && continue

        ix_na = tr._neighbour_a
        if (ix_na > ix) | isexternal(t._trigs[ix_na])
            push!(result, DelaunayEdge(getb(tr), getc(tr)))
        end
        ix_nb = tr._neighbour_b
        if (ix_nb > ix) | isexternal(t._trigs[ix_nb])
            push!(result, DelaunayEdge(geta(tr), getc(tr)))
        end
        ix_nc = tr._neighbour_c
        if (ix_nc > ix) | isexternal(t._trigs[ix_nc])
            push!(result, DelaunayEdge(geta(tr), getb(tr)))
        end
    end
    return result
end
