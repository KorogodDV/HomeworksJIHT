function solve_cubic(a, b, c, d)
    # convert to Ax³ + 3Bx² + 3Cx + D = 0
    A, B, C, D = promote(float(a), float(b), float(c), float(d)) ./ (1, 3, 3, 1)
    δ₁ = A * C - B * B
    δ₂ = A * D - B * C
    δ₃ = B * D - C * C
    d13 = δ₁ * δ₃
    d22 = δ₂ * δ₂
    Δ = 4 * d13 - d22
    nanvalue = zero(A) / zero(A)

    if Δ < 0
        At = Cb = Db = zero(A)  # A-tilde, C-bar, D-bar
        if B^3 * D ≥ A * C^3
            At, Cb, Db = A, δ₁, -2 * B * δ₁ + A * δ₂
        else
            At, Cb, Db = D, δ₃, -D * δ₂ + 2 * C * δ₃
        end
        T₀ = -copysign(At, Db) * sqrt(-Δ)
        T₁ = -Db + T₀
        p = cbrt(T₁ / 2)
        q = T₁ == T₀ ? -p : -Cb / p
        x₁ = Cb ≤ 0 ? p + q : -Db / (p^2 + q^2 + Cb)
        x, w = B^3 * D ≥ A * C^3 ? (x₁ - B, A) : (-D, x₁ + C)
        return (x/w, nanvalue, nanvalue)
    else
        δ₁ == δ₂ == δ₃ == 0 && return (-B/A, -B/A, -B/A)
        sΔ = sqrt(Δ)
        θA, θD = abs.(atan.((A*sΔ, 2*B*δ₁ - A*δ₂, D*sΔ, D*δ₂ - 2*C*δ₃)) ./ 3)
        sCA, sCD = sqrt.(.-min.((δ₁, δ₃), 0))
        x₁A, x₁D = 2 .* (sCA, sCD) .* cos.((θA, θD))
        x₃A, x₃D = .-((sCA, sCD)) .* (cos.((θA, θD)) .+ sqrt(3) .* sin.((θA, θD)))
        xlt = (x₁A + x₃A > 2 * B) ? x₁A : x₃A
        xst = (x₁D + x₃D < 2 * C) ? x₁D : x₃D
        xl, wl = xlt - B, A
        xs, ws = -D, xst + C
        E = wl * ws
        F = -xl * ws - wl * xs
        G = xl * xs
        xm, wm = C * F - B * G, C * E - B * F
        return (xs/ws, xm/wm, xl/wl)
    end
end