## pav(y) returns the isotonic regression of the predictor y.

function pav{T<:Real}(y::Vector{T})
    yy = similar(y, Float64)
    n = length(y)
    i = similar(y, Int)         # index
    l = similar(y, Int)         # length
    ci = 1
    i[ci] = l[ci] = 1
    yy[ci] = y[1]
    for j=2:n
        ci += 1
        i[ci] = j
        l[ci] = 1
        yy[ci] = y[j]
        while ci ≥ 2 && yy[max(ci-1,1)] ≥ yy[ci]
            nl = l[ci-1] + l[ci]
            yy[ci-1] += l[ci] * (yy[ci] - yy[ci-1]) / nl
            ci -= 1
            l[ci] = nl
        end
    end
    while n > 0
        for j=i[ci]:n
            yy[j] = yy[ci]
        end
        n = i[ci] - 1
        ci -= 1
    end
    yy
end
