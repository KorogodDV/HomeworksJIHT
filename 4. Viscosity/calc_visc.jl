include("./Viscosity/src/Viscosity.jl")

using .Viscosity
using Printf


calc_visc(gas::Gas, T_min::Real, T_max::Real) = reshape([T_min:T_max; visc.(Ref(gas), T_min:T_max)], :, 2)

CO₂ = Gas(44.009e-3, 3.996, 190)
CH₄ = Gas(16.043e-3, 3.822, 137)
O₂ = Gas(31.999e-3, 3.433, 113)


open("./Output data/O_2.txt", "w") do file
    for T in range(100, 1000)
        curr_visc = visc(O₂, T)
        write(file, "$T $curr_visc\n")
    end
end