"""
    rectangulardomain(Zb, Ze, Δr)

Generate initial mesh node coordinates for a rectangular domain ∈ {`Zb` to `Ze`} with
initial mesh step `Δr`.

See also: `rect_dom.m`
"""
function rectangulardomain(Zb::Complex, Ze::Complex, Δr)
    X = real(Ze) - real(Zb)
    Y = imag(Ze) - imag(Zb)

    # BUG?: Why are n and m in real and imag (x and y) different?
    n = ceil(Int, Y/Δr + 1)
    dy = Y/(n-1)
    m = ceil(Int, X/sqrt(Δr^2 - dy^2/4) + 1)
    dx = X/(m-1)

    vx = range(real(Zb), length=m, stop=real(Ze))
    vy = range(imag(Zb), length=n, stop=imag(Ze))

    tmp = ones(n)
    tmp[n] = 0.0

    # TODO: This looks _so_ matlaby
    y = repeat(vy, m)
    y .+= 0.5*dy*kron((1 .+ (-1).^(1:m))/2, tmp)
    y = [y; fill(imag(Zb), size(2:2:m))]
    x = reshape(repeat(vx', n), m*n)
    x = [x; ((2:2:m) .- 1)*dx .+ real(Zb)]

    # NOTE: Matlab values of `x` differ slightly because of Matlab's float handling and
    # transpose operator.
    # `sum(x)` is much closer between Julia and Matlab if the above line for `x` is:
    # x = [x' ((2:2:m) .- 1)'*dx .+ real(Zb)]'
    # x = reshape(x, length(x))

    return complex.(x, y)
end

"""
    diskdomain(R, Δr)

Generate initial mesh coordinates for a circular disk domain of radius `R` and center (0, 0)
for ``|z| < R``.

See also: `disk_dom.m`
"""
function diskdomain(R, Δr)
    h = Δr*sqrt(3)/2
    n = 1 + round(Int, R/h)
    Rn = (1:n)*R/n
    newnodes = [complex(zero(R))]
    f₀ = 0.0
    np = 6
    for ii = 1:n
        f = f₀ .+ range(0, stop=2π, length=np+1)
        xyn = Rn[ii]*cis.(f[1:end-1])
        append!(newnodes, xyn)
        f₀ += π/6/n
        np += 6
    end
    return newnodes
end
