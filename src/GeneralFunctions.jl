function distance(p1::Complex, p2::Complex)
    return sqrt((real(p2) - real(p1))^2 + (imag(p2) - imag(p1))^2)
end
function distance(p1::AbstractPoint2D, p2::AbstractPoint2D)
    return sqrt((getx(p2)-getx(p1))^2 + (gety(p2)-gety(p1))^2)
end
distance(e::DelaunayEdge) = distance(getb(e), geta(e))
function distance(p1arr::AbstractArray{IndexablePoint2D}, p2::IndexablePoint2D)
    return sqrt.((getx.(p1arr).-getx(p2)).^2 .+ (gety.(p1arr).-gety(p2)).^2)
end

"""
    longedge(edge, tolerance, geom2fcn)

Return true if `edge` has length greater than `tolerance`.
"""
function longedge(edge::DelaunayEdge, tolerance, geom2fcn::Function)
    return distance(geom2fcn(edge)...) > tolerance  # TODO: splat performance?
end
