function GaussianElimination!(x::Array{Float64, 2}, y::Array{Float64, 1})
    # define variables
    m, n = size(x)  # m: rows, n: cols
    δ = 1e-10

    if m != size(y)[1] || m != n
        exit(1)
    end

    for i in 1:(m-1)
        if abs(x[i, i]) < δ
            if (index = findNonZero(x[i:m, i])+i) == i
                continue
            else
                x[i, :], x[index, :] = x[index, :], x[i, :]
                y[i], y[index] = y[index], y[i]
            end
        end

        for j in (i+1):m
            y[j] -= (y[i] * x[j, i] / x[i, i])
            x[j, i:m] -= (x[i, i:m] .* x[j, i] ./ x[i, i])
        end
    end

    for i in m:(-1):2
        if abs(x[i, i]) < δ
            continue
        end

        for j in (i-1):(-1):1
            y[j] -= (y[i] * x[j, i] / x[i, i])
            x[j, i] = 0.0
        end
    end
end

function findNonZero(x::Array{Float64, 1})
    m = size(x)[1]
    δ = 1e-10

    for i in 1:m
        if abs(x[i]) < δ
            return i
        end
    end

    return 0
end
