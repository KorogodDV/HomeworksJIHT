function newtonsys(f, x, J; maxiter=50, xtol=1e-6, ftol=1e-6)
    x = float(copy(x))
    δx, y = similar.((x, x))
    for i in 1:maxiter
        y .= f(x)
        δx .= .- (J(x) \ y)
        x .+= δx

        norm(δx) < xtol && return x
        norm(y) < ftol && return x
    end
    error("Превышено число итераций.")
end