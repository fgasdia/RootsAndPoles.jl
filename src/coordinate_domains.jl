"""
    rectangulardomain(zbl, ztr, Δr)

Generate initial mesh node coordinates for a rectangular domain extending from the complex
bottom left corner `zbl` to the top right corner `ztr` with initial mesh step `Δr`.
"""
function rectangulardomain(zbl::Complex, ztr::Complex, Δr)
    rzbl, izbl = reim(zbl)
    rztr, iztr = reim(ztr)

    X = rztr - rzbl
    Y = iztr - izbl

    n = ceil(Int, Y/Δr)
    dy = Y/n

    ## dx = sqrt(Δr² - (dy/2)²), solved for equilateral triangle
    m = ceil(Int, X/sqrt(Δr^2 - dy^2/4))
    dx = X/m
    half_dx = dx/2

    T = promote_type(ComplexF64, typeof(zbl), typeof(ztr), typeof(Δr))
    v = Vector{T}()
    sizehint!(v, (m+1)*(n+1))

    shift = false  # we will displace every other line by dx/2
    for j = 0:n
        y = izbl + dy*j

        for i = 0:m
            x = rzbl + dx*i

            if shift && i == 0
                continue  # otherwise, we shift out of left bound
            elseif shift
                x -= half_dx
            end

            if i == m
                shift = !shift
            end

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
