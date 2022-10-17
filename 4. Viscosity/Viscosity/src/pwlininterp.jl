function hatfunc(t, k)
    n = size(t, 1)
    return function (x)
        if k ≥ 2 && t[k-1] ≤ x ≤ t[k]
            return (x - t[k-1]) / (t[k] - t[k-1])
        elseif k ≤ n-1 && t[k] ≤ x ≤ t[k+1]
            return (t[k+1] - x) / (t[k+1] - t[k])
        else  # x вне [t[1], t[end]] или неподходящий k
            return zero(x)
        end
    end
end

function pwlininterp(t, y)
    basis = [hatfunc(t, k) for k in 1:size(t, 1)]
    return x -> sum(y[k]*basis[k](x) for k in 1:size(y, 1))
end