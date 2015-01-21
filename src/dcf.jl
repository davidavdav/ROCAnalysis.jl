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
## For a Bayes Error Rate, scores scholf be log-likelihood ratios, and the Bayes decisions
## are made at threshold -plo.
## However, because this routines is also used for dcf(), we allow a separate setting of
## the threshold.
function ber{T1<:Real,T2<:Real}(tar::Vector{T1}, non::Vector{T1}, plo::T2, thres::T2=-plo)
    nmiss = nfa = 0
    for x in tar
        nmiss += x ≤ thres
    end
    for x in non
        nfa += x ≥ thres
    end
    ## the effective prior is sigmoid(plo), use sigmoid(-plo) for p_non for numerical stability
    return sigmoid(plo) * nmiss / length(tar) + sigmoid(-plo) * nfa / length(non)
end
## Both the prior log odds and the threshold may be an array.
## If the plo's are an array, the thresholds---if given---must be an array of the same size,
## because it doesn't really make sense to evaluate different cost functions with the same threshold. 
## For larger arrays of thres or plo, it is faster to use the Roc version below
ber{T1<:Real,T2<:Real}(tar::Vector{T1}, non::Vector{T1}, plo::Real, thres::Vector{T2}) = [ber(tar, non, plo, x) for x in thres]
function ber{T1<:Real,T2<:Real}(tar::Vector{T1}, non::Vector{T1}, plo::Vector{T2}, thres::Vector{T2}=-plo)
    size(plo) == size(thres) || error("Inconsistent vector lengths")
    [ber(tar, non, x...) for x in zip(plo,thres)]
end

## Use an _unclollapsed_ ROC for determining the Bayes error rate.
## r = roc(tar, non, collapse=false)

## We can't accurately compute this from a collapsed ROC, because
## the actual threshold might be somewhere among the set of collapsed scores. 

## Bayes false alarm and miss rates, we need this for the NBE plot.
## First define a helper function, scalar and array versions
function ber_famiss(r::Roc, field::Symbol, plo::Real, thres::Real=-plo)
    i = binsearch(thres, getfield(r, field)) + 1
    return sigmoid(-plo) * r.pfa[i], sigmoid(plo) * r.pmiss[i]
end
function ber_famiss{T<:Real}(r::Roc, f::Symbol, plo::Real, thres::Array{T})
    tuples = [ber_famiss(r, f, x) for x in thres]
    return [ [ x[y] for x in tuples] for y in 1:2]
end
function ber_famiss{T<:Real}(r::Roc, f::Symbol, plo::Array{T}, thres::Array{T}=-plo)
    size(plo) == size(plo) || error("Inconsistent vector lengths")
    tuples = [ber_famiss(r, f, x) for x in plo]
    return [ [ x[y] for x in tuples] for y in 1:2]
end

## actual errors
ber_famiss(r::Roc, plo::ArrayOrReal, thres=-plo) = ber_famiss(r, :θ, plo, thres)

ber(r::Roc, plo::Real, thres::Real=-plo) = +(ber_famiss(r, plo, thres)...)
ber{T<:Real}(r::Roc, plo::Real, thres::Array{T}) = [ber(r, plo, x) for x in thres]
ber{T<:Real}(r::Roc, plo::Array{T}, thres::Array{T}=-plo) = [ber(r, x...) for x in zip(plo,thres)]

## minimum ber is best computed using a Roc structure
minber_famiss(r::Roc, plo::ArrayOrReal) = ber_famiss(r, :llr, plo)
minber(r::Roc, plo::Real) = +(minber_famiss(r, plo)...)
minber{T<:Real}(r::Roc, plo::Array{T}) = [minber(r, x) for x in plo]

## compute ROC if scores are given
function minber{T<:Real}(tar::Vector{T}, non::Vector{T}, plo::ArrayOrReal)
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

## Allow a default DCF to be set
Base.copy(d) = DCF(d.ptar, d.cfa, d.cmiss)
function setdcf(;ptar=0.5, cfa=1, cmiss=1, d=DCF(ptar, cfa, cmiss))
    global global_dcf = copy(d)
end
## set a default DCF
setdcf()
getdcf() = global_dcf


dcf(tar::Vector, non::Vector; d::DCF=getdcf(), thres=-plo(d), norm=false) = applyfactor(d, ber(tar, non, plo(d), thres), norm)
dcf(tnt::TNT; d::DCF=getdcf(), thres=-plo(d), norm=false) = dcf(tnt.tar, tnt.non, d=d, thres=thres, norm=norm)
dcf(r::Roc; d::DCF=getdcf(), thres=-plo(d), norm=false) = applyfactor(d, ber(r, plo(d), thres), norm)

mindcf(r::Roc; d::DCF=getdcf(), norm=false) = applyfactor(d, minber(r, plo(d)), norm)
mindcf(tar::Vector, non::Vector; d::DCF=getdcf(), norm=false) = applyfactor(d, minber(roc(tar, non), plo(d)), norm)
mindcf(tnt::TNT; d::DCF=getdcf(), norm=false) = mindcf(tnt.tar, tnt.non, d=d, norm=norm)
