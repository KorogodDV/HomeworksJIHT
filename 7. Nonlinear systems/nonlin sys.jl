using Plots
using ForwardDiff
using LinearAlgebra
using LaTeXStrings

function broydensys(f, x, J; maxiter=1000, xtol=1e-10, ftol=1e-10)
    δx = float(similar(x))
    yp, yn = similar.((δx, δx))
    x = float(copy(x))
    B = float(copy(J))
    yn .= f(x)
    for i in 1:maxiter
        yp .= yn
        δx .= .- (B \ yp)
        x .+= δx
        yn .= f(x)

        norm(δx) < xtol && return x
        norm(yn) < ftol && return x
        
        g = B * δx
        B .+= (1 / dot(δx, δx)) .* (yn .- yp .- g) .* δx'
    end
    error("Превышено число итераций.")
end

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

create_P(T::Float64) = (V) -> 8T/(3V - 1)-3/V^2
create_μ(T::Float64) = (V) -> -T*log((3V-1)/(2ℯ^(-1/2))) + T/(3V - 1) - 9/(4V)

create_f(T::Float64) = V -> [create_P(T)(V[1]) - create_P(T)(V[2]), create_μ(T)(V[1]) - create_μ(T)(V[2])]

# Выбрал метод Ньютона, потому что метод Бройдена дает
# неправильный результат при T = 0.99.
function calc_equilibrium(T)
    f = create_f(T)
    J(x) = ForwardDiff.jacobian(f, x)
    x0 = [0.6, 3.0]
    # return broydensys(f, x0, J(x0))
    return newtonsys(f, x0, J)
end

T_range = 0.85:0.02:0.99
V_range = range(0.35, 5, 200)
V_eq = Vector{Float64}()
P_eq = Vector{Float64}()
plt = plot(; xlim=(0.3, 3.3), ylim=(0, 1.33), xlabel=L"V/V_c", ylabel=L"P/P_c")

for T in T_range
    V_eq_T = calc_equilibrium(T)
    push!(V_eq, V_eq_T...)
    push!(P_eq, create_P(T).(V_eq_T)...)
    plot!(plt, V_range, create_P(T).(V_range); label=L"T_r = %$T",
    line = (:dash, 2))
end

V_eq = reshape(V_eq, (2, :))
P_eq = reshape(P_eq, (2, :))
println(V_eq, P_eq)
plot!(plt, [V_eq[1, :], V_eq[2, :]], [P_eq[1, :], P_eq[2, :]]; 
label = [L"V_L" L"V_G"],
line = (4), marker = ([:circle :circle]))
display(plt)