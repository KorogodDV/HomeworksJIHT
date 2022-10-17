include("./Viscosity/src/Viscosity.jl")

using .Viscosity
using Printf


calc_visc(gas::Gas, T_min::Real, T_max::Real) = reshape([T_min:T_max; visc.(Ref(gas), T_min:T_max)], :, 2)

CO₂ = Gas(44.009 * 1e-3, 3.996, 190)
CH₄ = Gas(16.043 * 1e-3, 3.822, 137)
O₂ = Gas(31.999 * 1e-3, 3.433, 113)


file = open("D:/Homeworks/Теплофизика/HomeworksJIHT/3. Viscosity/Output data/O_2.txt", "w")
for T in range(100, 1000)
    curr_visc = visc(O₂, T)
    write(file, "$T $curr_visc\n")
end
close(file)