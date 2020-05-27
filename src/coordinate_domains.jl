"""
    rectangulardomain(Zb, Ze, Δr)

Generate initial mesh node coordinates for a rectangular domain ∈ {`Zb`, `Ze`} with
initial mesh step `Δr`.

`Ze` should be greater than `Zb`.
"""
function rectangulardomain(Zb::Complex, Ze::Complex, Δr)
    rZb, iZb = reim(Zb)
    rZe, iZe = reim(Ze)

    # todo: better error
    rZb > rZe && throw(ArgumentError("real part of Zb greater than Ze"))
    iZb > iZe && throw(ArgumentError("imag part of Zb greater than Ze"))

    X = rZe - rZb
    Y = iZe - iZb

    n = ceil(Int, Y/Δr + 1)
    dy = Y/(n-1)
    half_dy = dy/2

    m = ceil(Int, X/sqrt(Δr^2 - dy^2/4) + 1)
    dx = X/(m-1)

    # precalculate
    mn = m*n

    vlength = mn + div(m, 2)
    v = Vector{promote_type(typeof(Zb), typeof(Ze), typeof(Δr))}(undef, vlength)

    I = LinearIndices((n,m))
    on = false
    @inbounds for j in 0:m-1, i in 0:n-1
        # dx/j and dy/i flipped from what would seem "natural" to match Matlab meshgrid
        x = rZb + dx*j
        y = iZb + dy*i

        # `n` zeros then `n-1` ones and one zero, then repeat: `n` zeros, etc
        if (i+1) % n == 0
            on = !on
        elseif on && ((i+1) % n != 0)
            y += half_dy
        end
        v[I[i+1,j+1]] = complex(x, y)
    end

    idx = mn
    @inbounds for i in 1:2:(m-1)
        idx += 1

        # The "extra" rows or columns still cover the whole range from `Zb` to `Ze`
        tx = rZb + dx*i
        ty = iZb
        v[idx] = complex(tx, ty)
    end

    return v
end

"""
    diskdomain(R, Δr)

Generate initial mesh coordinates for a circular disk domain of radius `R` and
center (0, 0) for ``|z| < R`` and initial mesh step `Δr`.
"""
function diskdomain(R, Δr)
    h = Δr*sqrt(3)/2
    n = 1 + round(Int, R/h)
    R_n = R/n

    f₀ = 0.0
    np = 6

    f₀step = π/(6n)

    newnodes = Vector{complex(promote_type(typeof(R), typeof(Δr)))}(undef, 6*sum(1:n)+1)
    newnodes[1] = 0

    idx = 2
    for ii = 1:n
        # precalculate
        iiR_n = ii*R_n
        step = 2π/np

        for jj = 0:np-1
            f = f₀ + step*jj

            newnodes[idx] = iiR_n*cis(f)
            idx += 1
        end

        f₀ += f₀step
        np += 6
    end
    return newnodes
end
