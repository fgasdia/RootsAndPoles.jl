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
    triangulardomain(Za, Zb, Zc, Δr)

Generate initial mesh node coordinates for a *grid-aligned right* triangle domain
∈ {`Za`, `Zb`, `Zc`} with initial mesh step `Δr`.

This function generates a mesh grid for particular right triangles where two
sides of the triangle are aligned to the underlying real/imaginary grid. Examples
of such triangles are:

a ----- b
|     /
|   /
| /
c

where
    - a, b have greatest extent in x
    - a, c have greatest extent in y
"""
function triangulardomain(Za::Complex, Zb::Complex, Zc::Complex, Δr)
    rZa, iZa = reim(Za)
    rZb, iZb = reim(Zb)
    rZc, iZc = reim(Zc)

    #==
    # Check if this is a right triangle
    validtriangle = true
    if rZa == rZb == rZc
        validtriangle = false
    elseif iZa == iZb == iZc
        validtriangle = false
    elseif rZa == rZb
        iZa == iZc || iZb == iZc || (validtriangle = false)
    elseif rZa == rZc
        iZa == iZb || iZb == iZc || (validtriangle = false)
    elseif rZb == rZc
        iZa == iZb || iZa == iZc || (validtriangle = false)
    else
        validtriangle = false
    end
    validtriangle || throw(ArgumentError("`Za`, `Zb`, `Zc` do not define a grid-aligned right triangle"))

    iZa == iZb || ((Zb, Zc) = (Zc, Zb))
    rZb > rZa || ((Za, Zb) = (Zb, Za))
    ==#

    # Determine `dx` and `dy`
    X = rZb - rZa
    Y = abs(iZa - iZc)

    n = ceil(Int, Y/Δr + 1)
    dy = Y/(n-1)
    half_dy = dy/2

    m = ceil(Int, X/sqrt(Δr^2 - dy^2/4) + 1)
    dx = X/(m-1)

    # precalculate
    mn = m*n
    slope = Y/X

    vlength = mn + div(m, 2)  # BUG?
    T = promote_type(typeof(Za), typeof(Zb), typeof(Zc), typeof(Δr), Float64)
    v = Vector{T}()

    on = false
    for j = 0:m-1  # col
        for i = 0:n-1  # row (we're traversing down column)

            x = rZa + dx*j
            y = iZc + dy*i

            if (i+1) == n
                on = !on
            end
            if on
                y -= half_dy
            end

            if y >= (iZc + slope*(x - rZa))
                push!(v, complex(x, y))
            end
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
