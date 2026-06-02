"""
    union_unoriented_edges!(u, E...)

Empty and fill `u` with each unique unoriented edge of collections `E`. 

!!! warning

    No explicit type assert is performed for `u`.
"""
function union_unoriented_edges!(u, E...)
    empty!(u)
    for e in E
        for x in e
            DT.contains_unoriented_edge(x, UE) || push!(u, x)
        end
    end
end

function append_unoriented_edges!(E, edges)
    for e in edges
        DT.contains_unoriented_edge(e, E) || push!(E, e)
    end
end

"""
    dedupe!(v, tol=1e-9, digits=toldigits(tol))
    dedupe!(roots, poles; tol=1e-9, digits=toldigits(tol))

Modify `v` in place, retaining only unique values up to the specified number of `digits`
after the decimal place.
"""
dedupe!(v; tol=1e-9, digits=toldigits(tol)) = unique!(x -> round(x; digits) + zero(x), v)
function dedupe!(roots, poles; tol=1e-9, digits=toldigits(tol))
    dedupe!(roots; digits)
    dedupe!(poles; digits)
end

dedupe(v; tol=1e-9, digits=toldigits(tol)) = unique(x -> round(x; digits) + zero(x), v)

"""
    toldigits(tol)
    
Return 1 less digit after the decimal place than the tolerance `tol`.

For example, `toldigits(1e-9) == 8`.
"""
toldigits(tol) = -round(Int, log10(tol), RoundDown) - 1
