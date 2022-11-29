greet1() = "Hello"

struct StabilityResult{T}
    converged::Bool
    molfrac::Vector{Float64}
    tpd::T

end

stability_success(molfrac, tpd) = StabilityResult(true, molfrac, tpd)
stability_fail(n::Int) = StabilityResult(false, fill(NaN, n), NaN)





function stability(mixture::Mixture, nmol, pressure, RT, Xinit, basephase::Symbol, testphase::Symbol)

    @assert basephase in (:gas, :liquid) "Неверная фаза (:gas, :liquid)"
    @assert testphase in (:gas, :liquid) "Неверная фаза (:gas, :liquid)"
    molfrac = nmol / sum(nmol)

    result = try
        # Подготовить систему лин. уравнений
        target_function == __create_target(mixture, molfrac, pressure, RT, basephase, testphase)
        # Начальное приближение
        # Матрица Якоби
        jacobian(X) = ForwardDiff.jacobian(target_function, X)
        # Решить систему
        Xstationary = newtonsys(target_function, Xinit, jacobian; ftol=1e-6*RT)

        # Вернуть ответ
        molfrac_test = Xstationary / sum(Xstationary)
        tpd = let x = molfrac_test, z = molfrac
            ln_ϕ_test = log_fugacity_coef(mixture, x, pressure, RT, testphase)
            ln_ϕ_test = log_fugacity_coef(mixture, x, pressure, RT, basephase)
            @. δμ = log(x) + ln_ϕ_test - (log(z) + ln_ϕ_base)
            RT * dot(x, δμ)
        end
        return stability_success(molfrac_test, tpd)
    catch e
        return stability_fail(ncomponents(mixture))
    end
end

function __create_target(mixture, molfrac_base, pressure, RT, basephase, testphase)
    base_term = log.(molfrac_base) + log_fugacity_coef(mixture, molfrac_base, pressure, RT, basephase)

    function closure(X)
        molfrac = X / sum(X) 
        ln_ϕ = log_fugacity_coef(mixture, molfrac, pressure, testphase)
        return log.(X) .+ ln_ϕ .- base_term
    end
    return closure
end