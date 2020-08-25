"""
    rectangulardomain(Zb, Ze, Δr)

Generate initial mesh node coordinates for a rectangular domain ∈ {`Zb`, `Ze`} with
initial mesh step `Δr`.

`Ze` should be greater than `Zb`.
"""
function rectangulardomain(Zb::Complex, Ze::Complex, Δr)
    rZb, iZb = reim(Zb)
    rZe, iZe = reim(Ze)

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

    vlength = mn
    T = promote_type(typeof(Zb), typeof(Ze), typeof(Δr), Float64)
    v = Vector{T}()
    sizehint!(v, vlength)

    on = false
    @inbounds for j = 0:m-1, i = 0:n-1
        x = rZb + dx*j
        y = iZb + dy*i

        if (i+1) == n
            on = !on
        end
        if on
            y += half_dy
        end

        if y < iZe || (abs(y - iZe) < sqrt(eps(real(T))))
            # Can't just check y <= iZe because of floating point limitations
            push!(v, complex(x, y))
        end
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

    T = promote_type(typeof(R), typeof(Δr), Float64)
    newnodes = Vector{complex(T)}(undef, 6*sum(1:n)+1)
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
