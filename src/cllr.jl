## cllr.jl  Cost of the Log Likelihood Ratio
## (c) 2014--2015 David A. van Leeuwen
##
## Licensed under the MIT software license, see LICENSE.md

## Cllr for target and non-target scores
function cllr{T<:Real}(tar::Vector{T}, non::Vector{T})
    ct = cn = zero(T)
    for x in tar
        ct += softplus(-x)
    end
    for x in non
        cn += softplus(x)
    end
    (ct / length(tar) + cn / length(non)) / 2log(2)
end

## minimum Cllr: Cllr after optimal score-to-llr transformation
mincllr{T<:Real}(tar::Vector{T}, non::Vector{T}) = cllr(optllr(tar, non)...)
