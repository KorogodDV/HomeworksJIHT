using Plots

function intadapt(f, a, b, tol, xtol=eps(), fa=f(a), fb=f(b), m=(b-a)/2, fm=f(m))
    if a > b; a, b = b, a; end
    
    xl = (a + m)/2; fl = f(xl)  # расположение:
    xr = (m + b)/2; fr = f(xr)  # a -- xl -- m -- xr -- b
    
    T = Vector{Float64}(undef, 3)
    h = b - a
    T[1] = h * (fa + fb)/2
    T[2] = T[1]/2 + h/2 * fm
    T[3] = T[2]/2 + h/4 * (fl + fr)
    S = (4*T[2:end] - T[1:2]) / 3

    err = (S[2] - S[1]) / 15
    
    if abs(err) < tol * (1 + tol * abs(S[2]))
        Q = S[2]
        nodes = [a, xl, m, xr, b]
    else
        b - a ≤ xtol && error("Достигнут предел точности отрезка интегрирования `xtol`.")
        Ql, nodesl = intadapt(f, a, m, tol, xtol, fa, fm, xl, fl)
        Qr, nodesr = intadapt(f, m, b, tol, xtol, fm, fb, xr, fr)
        Q = Ql + Qr
        nodes = [nodesl; nodesr[2:end]]
    end
    return (Q, nodes)
end

function calc_1()
    a, b = - 1 // 3, 3
    f(x) = sin(100*x*exp(-x^2))
    Q, nodes = intadapt(f, a, b, 1e-6)
    return Q
end

function calc_2()
    a, b = 0, pi / 2
    f(x) = exp(x)*cos(x)
    res = Vector[]
    for i in -2:-1:-12
        Q, nodes = intadapt(f, a, b, 10.0^i)
        push!(res, [i, Q, length(nodes)])
    end
    return res
end

println("5.5.1:  ", calc_1())
ad_res = calc_2()
println(ad_res)
ex_res = (ℯ^(pi / 2) - 1) / 2

# x = 1:10; y = rand(10, 2) # 2 columns means two lines
# p = plot(x, y)
# z = rand(10)
# plot!(p, x, z)

# a, b = 0, 1.85
# foo(x) = x * sin(2x/(x-2))
# Q, nodes = intadapt(foo, a, b, 1e-4)
# plot(; layout=(2,1), xlabel=L"x", leg=:topleft, linewidth=2)
# plot!(foo, a, b; label="функция", linewidth=2, linecolor=:red)
# plot!(nodes, foo.(nodes);
#     label="узлы интегрирования",
#     seriestype=:sticks,
#     marker=(:o, 2, :lightblue),
#     linecolor=:blue,
# )
# histogram!(nodes;
#     label="распределение числа узлов",
#     subplot=2,
#     bins=range(a, b; length=9),
#     link=:x,
# )