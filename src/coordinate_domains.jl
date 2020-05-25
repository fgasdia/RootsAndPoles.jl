"""
    rectangulardomain(Zb, Ze, Δr)

Generate initial mesh node coordinates for a rectangular domain ∈ {`Zb`, `Ze`} with
initial mesh step `Δr`.
"""
function rectangulardomain(Zb::Complex, Ze::Complex, Δr)
    rZb, iZb = reim(Zb)
    rZe, iZe = reim(Ze)

    X = rZe - rZb
    Y = iZe - iZb

    n = ceil(Int, Y/Δr + 1)
    dy = Y/(n-1)
    m = ceil(Int, X/sqrt(Δr^2 - dy^2/4) + 1)
    dx = X/(m-1)

    vx = range(rZb, stop=rZe, length=m)
    vy = range(iZb, stop=iZe, length=n)

    y = repeat(vy, m)

    half_dy = dy/2
    on = false
    for i in eachindex(y)
        if i % n == 0
            on = !on
        elseif on
            y[i] += half_dy
        end
    end

    y = vcat(y, fill(iZb, div(m,2)))  # TODO: why?

    x = repeat(vx, inner=n)
    x = vcat(x, (1:2:(m-1))*dx .+ rZb)

    return complex.(x, y)
end

"""
    diskdomain(R, Δr)

Generate initial mesh coordinates for a circular disk domain of radius `R` and center (0, 0)
for ``|z| < R``.
"""
function diskdomain_new(R, Δr)
    h = Δr*sqrt(3)/2
    n = 1 + round(Int, R/h)
    R_n = R/n

    f₀ = 0.0
    np = 6

    f₀step = π/(6n)

    newnodes = Vector{complex(typeof(R))}(undef, 6*sum(1:n)+1)
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
