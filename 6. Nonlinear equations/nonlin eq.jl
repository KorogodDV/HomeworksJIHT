using Plots

function itproot(f, x₁, x₂; xtol=eps(), ftol=eps(), κ₁=0.1, κ₂=2, n₀=1)
    if x₁ > x₂; x₁, x₂ = x₂, x₁; end
    y₁, y₂ = f(x₁), f(x₂)
    y₁ * y₂ > 0 && error("Функция должна иметь разные знаки в концах отрезка")
    y₁ == 0 && return x₁
    y₂ == 0 && return x₂
    
    nbisect = ceil(Int, log2((x₂-x₁)/xtol))
    maxiter = nbisect + n₀
    brackorig = x₂ - x₁

    for i in 1:maxiter
        # interpolate
        xf = (y₂*x₁ - y₁*x₂)/(y₂ - y₁)

        # truncate
        xmid = (x₁ + x₂)/2
        σ = sign(xmid - xf)
        δ = κ₁ * (x₂ - x₁)^κ₂ / brackorig
        xt = δ ≤ abs(xmid - xf) ? xf + copysign(δ, σ) : xmid
        
        # project
        r = xtol * 2.0^(maxiter - i) - (x₂ - x₁)/2
        xnew = abs(xt - xmid) ≤ r ? xt : xmid - copysign(r, σ)

        ynew = f(xnew)
        if sign(y₂) == sign(ynew)
            x₂, y₂ = xnew, ynew
        elseif sign(y₁) == sign(ynew)
            x₁, y₁ = xnew, ynew
        else  # ynew == 0
            return xnew
        end

        abs(ynew) < ftol && (x₁ + x₂)/2
        abs(x₂ - x₁) < xtol && (x₁ + x₂)/2
    end
    return (x₁ + x₂)/2
end

function rachford_rice_solve(z, K; xtol = 2e-10)
    sum(z) != 1 && error("Sum of the z elements must be equal to 1")
    F(G) = sum(z .* (K .- 1) ./ (G * (K .- 1) .+ 1))
    return itproot(F, 1 / (1 - maximum(K)), 1 / (1 - minimum(K)); xtol = xtol, n₀=1)
end

function solution_1()
    z = [0.9, 0.1]
    K = [1.5, 0.01]
    F(G) = sum(z .* (K .- 1) ./ (G * (K .- 1) .+ 1))
    G_root = rachford_rice_solve(z, K)
    # p = plot(F, G_root - 0.01, G_root + 0.01)
    p = plot(F, 0.6, 0.9; label = "F(G)")
    scatter!(p, [G_root], [F(G_root)]; label = "root")
    xlabel!("G")
    ylabel!("F(G)")
    return G_root, p
end

function solution_2()
    z = [0.2463, 0.2208, 0.2208, 0.3121]
    K = [40, 25, 0.6, 0.005]
    F(G) = sum(z .* (K .- 1) ./ (G * (K .- 1) .+ 1))
    G_root = rachford_rice_solve(z, K)
    # p = plot(F, G_root - 0.01, G_root + 0.01)
    p = plot(F, 0.3, 0.8; label = "F(G)")
    scatter!(p, [G_root], [F(G_root)]; label = "root")
    xlabel!("G")
    ylabel!("F(G)")
    return G_root, p
end 

root, p = solution_2()
println(root)
display(p)