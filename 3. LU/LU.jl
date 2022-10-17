using LinearAlgebra

function backwardsub!(x::Vector{Float64}, U::AbstractMatrix, b::AbstractVector)
    size(U, 1) != size(U, 2) && throw(ArgumentError("Matrix should be upper triangular and invertible"))
    !LinearAlgebra.istriu(U) && throw(ArgumentError("Matrix should be upper triangular and invertible"))
    any(≈(0), diag(U)) && throw(ArgumentError("Matrix should be invertible"))
    x[end] = b[end] / U[end]
    for i in size(U, 1)-1:-1:1
        x[i] = (b[i] - LinearAlgebra.dot(U[i, i+1:end], x[i+1:end])) / U[i, i]
    end
end

"This function solves the linear system Ux=b for upper triangular invertible REAL matrix U using backward substitution algorithm"
function backwardsub(U::AbstractMatrix, b::AbstractVector)
    x = zeros(size(b))
    backwardsub!(x, U, b)
    return x
end


function forwardsub!(x::Vector{Float64}, L::AbstractMatrix, b::AbstractVector)
    size(L, 1) != size(L, 2) && throw(ArgumentError("Matrix should be lower triangular and invertible"))
    !LinearAlgebra.istril(L) && throw(ArgumentError("Matrix should be lower triangular and invertible"))
    any(≈(0), diag(L)) && throw(ArgumentError("Matrix should be invertible"))
    x[1] = b[1] / L[1]
    for i in 2:size(L, 1)
        x[i] = (b[i] - LinearAlgebra.dot(L[i, 1:i-1], x[1:i-1])) / L[i, i]

    end
end

"This function solves the linear system Ux=b for lower triangular REAL invertible matrix U using forward substitution algorithm"
function forwardsub(L::AbstractMatrix, b::AbstractVector)
    x = zeros(size(b))
    forwardsub!(x, L, b)
    return x
end


"В обеих tridiagsolve будут обозначения не такие, как в задании, потому что я просто скопировал и адаптировал уже написанный на вычматах питоновский
код, а там у нас были другие обозначения. Надеюсь, это не большая проблема."

function tridiagsolve(c::AbstractVector, b::AbstractVector, a::AbstractVector, d::AbstractVector)
    L = length(a)
    length(d) != L && throw(ArgumentError("All arguments must have the same length"))
    abs(b[1]) < abs(a[1]) + abs(c[1]) && throw(ArgumentError("Матрица должна обладать диагональным преобладанием"))
    alpha = zeros(L)
    beta = zeros(L)
    alpha[1] = -1 * a[1] / b[1]
    beta[1] = d[1] / b[1]
    u = zeros(L)
    for i in 2:L
        abs(b[i]) < abs(a[i]) + abs(c[i]) && throw(ArgumentError("Матрица должна обладать диагональным преобладанием"))
        alpha[i] = -1 * a[i] / (b[i] + c[i] * alpha[i - 1])
        beta[i] =  (d[i] - c[i] * beta[i - 1]) / (b[i] + c[i] * alpha[i - 1])
    end
    u[L] = (d[L] - c[L] * beta[L - 1]) / (b[L] + c[L] * alpha[L - 1])
    for i in L-1:-1:1
        u[i] = alpha[i] * u[i + 1] + beta[i]
    end
    return u
end

function tridiagsolve(A::Tridiagonal, d::AbstractVector)
    L = size(A, 1)
    a = [0; A[2:L+1:end - L]]
    b = A[1:L+1:end]
    c = [A[L+1:L+1:end - 1]; 0]
    return tridiagsolve(a, b, c, d)
end

function solution_a()
    U = [8 9 4 -1; 0 4 1 0; 0 0 -1 6; 0 0 0 11]
    b = [9, 3, -1, 2]
    x = backwardsub(U, b)
    return x, b - U * x
end

function solution_b()
    T = [-2 1 0 0 0; 1 -2 1 0 0; 0 1 -2 1 0; 0 0 1 -2 1; 0 0 0 1 -2]
    f = [1, 1, 1, 1, 1]
    x = tridiagsolve(Tridiagonal(T), f)
    return x, f - T * x
end

function solution_c()
    A = [1 8 -3 9; 0 4 10 -2; 8 2 -5 1; 3 1 6 2]
    b = [3, 6, 1, 4]
    L, U, p = lu(A)
    x = backwardsub(U, forwardsub(L, b[p]))
    return x, b - A * x
end

# U = [8 9 4 -1; 0 4 1 0; 0 0 -1 6; 0 0 0 11]
# L = copy(U')
# T = [-2 1 0 0 0; 1 -2 1 0 0; 0 1 -2 1 0; 0 0 1 -2 1; 0 0 0 1 -2]
# b = [9, 3, -1, 2]
# f = [1, 1, 1, 1, 1]
# a = similar(b, Float64)
# println(backwardsub(U, b))
# println(tridiagsolve([0, 1, 1, 1, 1], [-2, -2, -2, -2, -2], [1, 1, 1, 1, 0], f))
# println(tridiagsolve(Tridiagonal(T), f))
# println(T \ f)
println(solution_b())
