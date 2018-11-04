"""
Distance from Complex `p1` to `p2`.
"""
function distance(p1::Complex, p2::Complex)
    sqrt((real(p2) - real(p1))^2 + (imag(p2) - imag(p1))^2)
end
distance(p1::AbstractPoint2D, p2::AbstractPoint2D) = sqrt((getx(p2)-getx(p1))^2 + (gety(p2)-gety(p1))^2)
distance(e::DelaunayEdge) = distance(getb(e), geta(e))
distance(p1arr::AbstractArray{IndexablePoint2D}, p2::IndexablePoint2D) = sqrt.((getx.(p1arr).-getx(p2)).^2 + (gety.(p1arr).-gety(p2)).^2)

"""
Return true if `edge` has length greater than `tolerance`.
"""
function longedge(edge::DelaunayEdge, tolerance, geom2fcn)
    distance(geom2fcn(edge)...) > tolerance
end
