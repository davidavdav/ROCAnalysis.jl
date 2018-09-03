## cllr.jl  Cost of the Log Likelihood Ratio
## (c) 2014--2015 David A. van Leeuwen
##
## Licensed under the MIT software license, see LICENSE.md

## Cllr for target and non-target scores
"""
`cllr(tar, non)` computes the cost of the log-likelihood-ratio, for log-likelihood-ratio (llr) scores
`tar` and `non` (target scores and non-target scores).  Target an non-target llr's are those for which 
the numerator `H1` and denominator `H2` hypothesis are actually true in the log-likelihood-ratio 

    llr(x) = p(x|H1) / p(x|H2)

where `x` is the test data for a system trying to disciminate between `H1` and `H2` given `x`.  

The Cllr of a perfect system, assigning a llr of `+Inf` t otarget scores and `-Inf` to non-target scores
is 0, for a system that is indifferent, producing a `llr=0` for every input `x` the `Cllr=1`.  
Please note that, for badly calibrated systems, `Cllr>1`.  The units of Cllr are measured in _bits_, 
and cllr can be seen as the average amount of information per trial that is _not_ extracted from the 
data. 

Cllr measures the _calibration_ as well as _discrimination_ ability of a system.  Discrimination 
entails the ability to produces target scores that are, typically, much higher than non-target scores.  
Calibration entails that for every individual `x` an optimal Bayes decision can be made if a 
cost function is known.  
"""
function cllr(tar::Vector{T}, non::Vector{T}) where T<:Real
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
"""
`mincllr(tar, non)` computes the `cllr` of target and non-target scores after an optimal
transformation of the data.  This allows for measuring the disrimination performance of
a system in the units of `cllr`, which are bits. 
"""
mincllr(tar::Vector{T}, non::Vector{T}) where {T<:Real} = cllr(optllr(tar, non)...)
