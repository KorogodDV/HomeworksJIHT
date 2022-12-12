using CTSolverStability
using Plots

c1 = VanDerWaalsComponent(2.283 / 10, 0.04278 / 1000, 4.5992e6, 190.56 * CTSolverStability.GAS_CONSTANT_SI)
c2 = VanDerWaalsComponent(5.562 / 10, 0.0638 / 1000, 4.88e6, 305 * CTSolverStability.GAS_CONSTANT_SI)
c3 = VanDerWaalsComponent(8.779 / 10, 0.08445 / 1000, 4.26e6, 370 * CTSolverStability.GAS_CONSTANT_SI)
c4 = VanDerWaalsComponent(14.66 / 10, 0.1226 / 1000, 3.796e6, 425.1 * CTSolverStability.GAS_CONSTANT_SI)
c5 = VanDerWaalsComponent(19.26 / 10, 0.146 / 1000, 3.38e6, 470 * CTSolverStability.GAS_CONSTANT_SI)
c6 = VanDerWaalsComponent(24.71 / 10, 0.1735 / 1000, 3.03e6, 508 * CTSolverStability.GAS_CONSTANT_SI)
n2 = VanDerWaalsComponent(1.370 / 10, 0.0387 / 1000, 3.39e6, 126 * CTSolverStability.GAS_CONSTANT_SI)

mixture = VanDerWaalsMixture([c1, c2, c3, c4, c5, c6, n2])
molfrac = [0.9430, 0.0270, 0.0074, 0.0049, 0.0027, 0.010, 0.0140]

pressures = range(1e5, 75e5; length=100)
temperatures = range(150, 250; length=150)

data = Vector{Any}[]
open("CTSolverStability/examples/data.tsv", "w") do io
    for P in pressures, T in temperatures
        RT = T * CTSolverStability.GAS_CONSTANT_SI
        Kinit = CTSolverStability.michelson.(mixture.components, P, RT)
        results = [
            CTSolverStability.stability(mixture, molfrac, P, RT, molfrac ./ Kinit, :gas, :liquid)
            CTSolverStability.stability(mixture, molfrac, P, RT, molfrac .* Kinit, :liquid, :gas)
        ]
        
        push!(data, Any[T, P, Int(results[1].converged), results[1].tpd, Int(results[2].converged), results[2].tpd])
        println(io, join(Any[T, P, Int(results[1].converged), results[1].tpd, 
        Int(results[2].converged), results[2].tpd
        ], '\t'))
    end
end

err, stab, unstab = [[], []], [[], []], [[], []]
for i in eachindex(data)
    curr_data = data[i]
    if curr_data[3] * curr_data[5] != 1
        push!(err[1], Float64(curr_data[1]))
        push!(err[2], Float64(curr_data[2]))
    elseif maximum(curr_data[[4, 6]]) <= -1e-11
        push!(unstab[1], Float64(curr_data[1]))
        push!(unstab[2], Float64(curr_data[2]))
    else
        push!(stab[1], Float64(curr_data[1]))
        push!(stab[2], Float64(curr_data[2]))
    end
end

p = scatter(stab[1], stab[2]; label="stable", xlabel="T, K", ylabel="P, Pa", marker = (:green, 2, :rect))
scatter!(p, unstab[1], unstab[2]; label="unstable", marker = (:purple, 2, :rect))
scatter!(p, err[1], err[2]; label="error", marker = (:red, 2, :rect))
display(p)
savefig(p, "CTSolverStability/examples/PT.png")

# Форма получившейся кривой позволяет оценить критическую точку просто как "пиковую" (самую дальнюю от центра) 
# на этой кривой. 
println("Критическая точка: T = $(unstab[1][end]) K, P = $(unstab[2][end]) Pa")