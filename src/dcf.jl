## dcf.jl  decision cost function support for the ROC package.
## (c) 2015 David A. van Leeuwen
##
## Licensed under the MIT software license, see LICENSE.md
##

## effective prior odds
oeff(ptar=0.5; cfa=1, cmiss=1) = (cmiss ./ cfa) .* (ptar ./ (1-ptar))
oeff(dcf::DCF) = oeff(dcf.ptar; cfa=dcf.cfa, cmiss=dcf.cmiss)
## effective prior
peff(ptar=0.5; cfa=1, cmiss=1) = 1 ./ (1 + (cfa ./ cmiss) * ((1-ptar) ./ ptar))
peff(dcf::DCF) = peff(dcf.ptar, cfa=dcf.cfa, cmiss=dcf.cmiss)
## prior log odds
plo(ptar=0.5; cfa=1, cmiss=1) = logit(ptar) + log(cmiss ./ cfa)
plo(dcf::DCF) = plo(dcf.ptar; cfa=dcf.cfa, cmiss=dcf.cmiss)

## Decision Cost Function: dcf = p_tar * C_miss * p_miss + (1-p_tar) * C_fa * p_fa
## Compute this through the Bayes error rate, ber

## Basic definition of what Niko calls Bayes Error Rate for a single prior log-odds
## count errors _at_ the decision boundary both in fa and miss
function ber{T1<:Real,T2<:Real}(tar::Vector{T1}, non::Vector{T1}, plo::T2)
    nmiss = nfa = 0
    for x in tar
        nmiss += x ≤ -plo
    end
    for x in non
        nfa += x ≥ -plo
    end
    p = sigmoid(plo)            # effective prior
    return p * nmiss / length(tar) + (1-p) * nfa / length(non)
end
## For larger arrays plo, it is faster to use the Roc version beow
ber{T1<:Real,T2<:Real}(tar::Vector{T1}, non::Vector{T1}, plo::Vector{T2}) = [ber(tar, non, x) for x in plo]

## Use an _unclollapsed_ ROC for determining the Bayes error rate.
## r = roc(tar, non, collapse=false)

## We can't accurately compute this from a collapsed ROC, because
## the actual threshold might be somewhere among the set of collapsed scores. 
function ber(r::Roc, plo::Real)
    i = binsearch(-plo, r.θ) + 1
    p = sigmoid(plo)
    return (1-p) * r.pfa[i] + p * r.pmiss[i]
end
ber{T<:Real}(r::Roc, plo::Vector{T}) = [ber(r, x) for x in plo]

## minimum ber is best computed using a Roc structure
function minber(r::Roc, plo::Real)
    i = binsearch(-plo, r.llr) + 1
    p = sigmoid(plo)
    return (1-p) * r.pfa[i] + p * r.pmiss[i]
end
minber{T<:Real}(r::Roc, plo::Vector{T}) = [minber(r, x) for x in plo]
## compute ROC if scores are given
function minber{T<:Real}(tar::Vector{T}, non::Vector{T}, plo::Real)
    r = roc(tar, non)
    return minber(r, plo)
end
function minber{T1<:Real,T2<:Real}(tar::Vector{T1}, non::Vector{T1}, plo::Vector{T2})
    r = roc(tar, non)
    return minber(r, plo)
end

## factor to normalize the Bayes error rate, for normalized bayes error rate (or normalized cdet)
normfactor(plo) = 1 + exp(abs(plo))

## factor to convert Bayes error rates to actual / minimum costs.  
costfactor(d::DCF) = d.ptar .* d.cmiss + (1-d.ptar) .* d.cfa

## to norm or not to norm
function applyfactor(d::DCF, er, norm::Bool)
    lo = plo(d)
    if norm
        return normfactor(lo) .* er
    else
        return costfactor(d) .* er
    end
end

## Decision cost functions

dcf(tar::Vector, non::Vector, d::DCF; norm=false) = applyfactor(d, ber(tar, non, plo(d)), norm)
dcf(tnt::TNT, d::DCF; norm=false) = dcf(tnt.tar, tnt.non, d, norm=norm)
dcf(r::Roc, d::DCF; norm=false) = applyfactor(d, ber(r, plo(d)), norm)

mindcf(r::Roc, d::DCF; norm=false) = applyfactor(d, minber(r, plo(d)), norm)
mindcf(tar::Vector, non::Vector, d::DCF; norm=false) = applyfactor(d, minber(roc(tar, non), plo(d)), norm)
mindcf(tnt::TNT, d::DCF; norm=false) = mindcf(tnt.tar, tnt.non, d, norm=norm)
