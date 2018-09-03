## pav.jl  Pool Adjacent Violators algorithm for isotonic regression
## (c) 2015 David A. van Leeuwen
##
## Licensed under the MIT software license, see LICENSE.md
## This code is largely based on code from the Bosaris toolkit

## pav(y) returns the isotonic regression of the predictor y.
function pav(y::Vector{T}) where T<:Real
    yy = similar(y, Float64)
    n = length(y)
    i = similar(y, Int)         # index
    l = similar(y, Int)         # length
    ci = 1
    i[ci] = l[ci] = 1
    yy[ci] = y[1]
    for j in 2:n
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
        for j in i[ci]:n
            yy[j] = yy[ci]
        end
        n = i[ci] - 1
        ci -= 1
    end
    return yy
end

## return the optimal log-likelihood ratios for target and non-target scores
function optllr(tar::Vector{T}, non::Vector{T}; laplace=true) where T<:Real
    ntar = length(tar)
    nnon = length(non)
    ntar > 0 && nnon > 0 || error("Lenghts must be nonzero")
    o = sortperm([tar; non])
    pideal = zeros(ntar + nnon + 4laplace)
    ## set pideal=1 for target scores
    for i in 1:ntar+nnon
        if o[i] ≤ ntar         # tar sorted before non, a target score
            pideal[i+2laplace] = 1.0
        end
    end
    if laplace
        pideal[1] = pideal[end-1] = 1.0
    end
    ## optimal posterior, given monotonicity constraints
    popt = pav(pideal)          # optimal posterior
    postlo = logit.(popt)[1+2laplace:end-2laplace] # posterior log odds
    priorlo = log(length(tar) / length(non))
    llrs = Array{T}(undef, ntar + nnon)
    llrs[o] = postlo .- priorlo
    llrs[1:ntar], llrs[ntar+1:end]
end
