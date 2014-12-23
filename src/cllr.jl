## (c) 2014 David A. van Leeuwen

## Cllr for target and non-target scores
function cllr{T<:Real}(tar::Vector{T}, non::Vector{T})
    ct = cn = zero(T)
    for x in tar
        ct += softplus(-x)
    end
    for x in non
        cn += softplus(x)
    end
    (ct / length(tar) + cn / length(non)) / (2log(2))
end
