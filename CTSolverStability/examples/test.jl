using CTSolverStability

c1 = VanDerWaalsComponent(2.283 / 10, 0.04278 / 1000, 4.5992e6, 190.56 * CTSolverStability.GAS_CONSTANT_SI)
c4 = VanDerWaalsComponent(14.66 / 10, 0.1226 / 1000, 3.796e6, 425.1 * CTSolverStability.GAS_CONSTANT_SI)

mixture = VanDerWaalsMixture([c1, c4])
molfrac = [0.8, 0.2]

pressures = range(1e5, 100e5; length=5)
temperatures = range(150, 500; length=5)

io = stdout

for P in pressures, T in temperatures
    RT = T * CTSolverStability.GAS_CONSTANT_SI
    Kinit = CTSolverStability.michelson.(mixture.components, P, RT)
    results = [
        CTSolverStability.stability(mixture, molfrac, P, RT, molfrac ./ Kinit, :gas, :liquid)
        CTSolverStability.stability(mixture, molfrac, P, RT, molfrac .* Kinit, :liquid, :gas)
    ]

    println(io, join(Any[P, T, Int(results[1].converged), results[1].tpd, 
    Int(results[2].converged), results[2].tpd
    ], '\t'))
end