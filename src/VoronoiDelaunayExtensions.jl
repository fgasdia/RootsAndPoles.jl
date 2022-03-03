"""
    IndexablePoint2D <: AbstractPoint2D

A special 2D point type, compatible with `VoronoiDelaunay`, that carries an
identifying `index` field.

By default, `index` is -1.

See also: `setindex!`, `getindex`
"""
mutable struct IndexablePoint2D <: AbstractPoint2D
    x::Float64
    y::Float64
    index::Int
end
IndexablePoint2D(x, y) = IndexablePoint2D(x, y, -1)
getx(p::IndexablePoint2D) = p.x
gety(p::IndexablePoint2D) = p.y
Base.:+(p1::IndexablePoint2D, p2::IndexablePoint2D) = IndexablePoint2D(getx(p1)+getx(p2), gety(p1)+gety(p2), -1)
Base.:/(p::IndexablePoint2D, n) = IndexablePoint2D(getx(p)/n, gety(p)/n, -1)

"""
    getindex(p::IndexablePoint2D)

Return `index` field of `p`.
"""
Base.getindex(p::IndexablePoint2D) = p.index

"""
    setindex!(p::IndexablePoint2D)

Set the `index` field of `p`.

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
Base.setindex!(p::IndexablePoint2D, v::Int) = setfield!(p, :index, v)

"""
    same(e1::DelaunayEdge{IndexablePoint2D}, e2::DelaunayEdge{IndexablePoint2D})

Return `true` if edge `e1` is the same as `e2`, regardless of direction.
Otherwise, return `false`. This is similar to `==`.

To check if two edges are identical in the sense that no program could
distinguish them, use `===`.
"""
@inline function same(e1::DelaunayEdge{IndexablePoint2D},e2::DelaunayEdge{IndexablePoint2D})
    e1 === e2 && return true

    # if they're reversed
    e1a, e1b = geta(e1).index, getb(e1).index
    e2a, e2b = geta(e2).index, getb(e2).index
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
        @inbounds edgei = v[i]
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
