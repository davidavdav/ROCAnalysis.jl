## dcf.jl  decision cost function support for the ROC package.
## (c) 2015 David A. van Leeuwen
##
## Licensed under the MIT software license, see LICENSE.md
##

## effective prior odds
"""
`oeff()`, the effective prior odds given cost function parameters.  Arguments are either:

 - `ptar, cfa, cmiss`: Scalars or Vectors of the prior of a target, the cost of a false positive and the cost of a false negative

 - `dcf::DCF`: a `DCF` object containing `ptar`, `cfa` and `cmiss`.
"""
oeff(ptar=0.5; cfa=1, cmiss=1) = (cmiss ./ cfa) .* (ptar ./ (1-ptar))
oeff(dcf::DCF) = oeff(dcf.ptar; cfa=dcf.cfa, cmiss=dcf.cmiss)
## effective prior
"""
`peff()`, the effective prior given cost function parameters. Arguments are either:

 - `ptar, cfa, cmiss`: Scalars or Vectors of the prior of a target, the cost of a false positive and the cost of a false negative

 - `dcf::DCF`: a `DCF` object containing `ptar`, `cfa` and `cmiss`.
"""
peff(ptar=0.5; cfa=1, cmiss=1) = 1 ./ (1 + (cfa ./ cmiss) * ((1-ptar) ./ ptar))
peff(dcf::DCF) = peff(dcf.ptar, cfa=dcf.cfa, cmiss=dcf.cmiss)
## prior log odds
"""
`plo()`, the prior log-odds given cost function parameters. Arguments are either:

 - `ptar, cfa, cmiss`: Scalars or Vectors of the prior of a target, the cost of a false positive and the cost of a false negative

 - `dcf::DCF`: a `DCF` object containing `ptar`, `cfa` and `cmiss`.
"""
plo(ptar=0.5; cfa=1, cmiss=1) = logit.(ptar) .+ log.(cmiss ./ cfa)
plo(dcf::DCF) = plo(dcf.ptar; cfa=dcf.cfa, cmiss=dcf.cmiss)

## Decision Cost Function: dcf = p_tar * C_miss * p_miss + (1-p_tar) * C_fa * p_fa
## Compute this through the Bayes error rate, ber

## Basic definition of what Niko calls Bayes Error Rate for a single prior log-odds
## count errors _at_ the decision boundary both in fa and miss
## For a Bayes Error Rate, scores schould be log-likelihood ratios, and the Bayes decisions
## are made at threshold -plo.
## However, because this routines is also used for dcf(), we allow a separate setting of
## the threshold.
"""
`ber()`, the Bayes Error Rate.  This computes the expected error rate given a set of supervised trials,
the log-likelihood-ratio scores `tar` and `non`, and a threshold based on a given a cost function
specified by its prior log-odds (see `plo()`).

If the scores `tar` and `non` are well-calibrated log-likelihood-ratios, the optimum decision threshold
that minimizes the decision cost function `dcf` is surprisingly simple,

    θ = -plo(dcf).

This function computes the actual cost for such decisions, given the test scores.

Arguments:

 - `tar::Vector`, `non::Vector`, `plo`: target and non-target llr
scores, and prior log odds.  Here, `plo` may be a vector, resulting in
multiple Bayes error rates.

 - `r::Roc`, `plo`: a `Roc` structure computed with `roc(;collapse=false)`.
It is mandatory that the Roc object is not collapsed, because the
actual cost may occur for a threshold between two "corner" points of
the ROC curve.  A `Roc` structure is useful for computing the minimum
Bayes Error Rate for multiple cost functions---in this case the `Roc`
may be `collapse`d.

 - `; thres`: an optional threshold, overriding the theoretical
optimum described above.  This can be used if the scores are not
correctly calibrated.
"""
function ber(tar::Vector{T1}, non::Vector{T1}, plo::Real, thres::Real=-plo) where T1<:Real
    nmiss = nfa = 0
    for x in tar
        nmiss += x ≤ thres
    end
    for x in non
        nfa += x ≥ thres
    end
    ## the effective prior is sigmoid(plo), use sigmoid(-plo) for p_non for numerical stability
    return sigmoid.(plo) * nmiss / length(tar) + sigmoid.(-plo) * nfa / length(non)
end
## Both the prior log odds and the threshold may be an array.
## If the plo's are an array, the thresholds---if given---must be an array of the same size,
## because it doesn't really make sense to evaluate different cost functions with the same threshold.
## For larger arrays of thres or plo, it is faster to use the Roc version below
ber(tar::Vector{T1}, non::Vector{T1}, plo::Real, thres::Vector{T2}) where {T1<:Real,T2<:Real} = [ber(tar, non, plo, x) for x in thres]
function ber(tar::Vector{T1}, non::Vector{T1}, plo::Vector{T2}, thres::Vector{T2}=-plo) where {T1<:Real,T2<:Real}
    size(plo) == size(thres) || error("Inconsistent vector lengths")
    [ber(tar, non, x...) for x in zip(plo, thres)]
end

## Use an _unclollapsed_ ROC for determining the Bayes error rate.
## r = roc(tar, non, collapse=false)

## We can't accurately compute this from a collapsed ROC, because
## the actual threshold might be somewhere among the set of collapsed scores.

## Bayes false alarm and miss rates, we need this for the NBE plot.
## First define a helper function, scalar and array versions
function ber_famiss(r::Roc, field::Symbol, plo::Real, thres::Real=-plo)
    i = binsearch(thres, getfield(r, field)) + 1
    return sigmoid.(-plo) * r.pfa[i], sigmoid.(plo) * r.pmiss[i]
end
function ber_famiss(r::Roc, f::Symbol, plo::Real, thres::Array{T}) where T<:Real
    tuples = [ber_famiss(r, f, x) for x in thres]
    return [ [ x[y] for x in tuples] for y in 1:2]
end
function ber_famiss(r::Roc, f::Symbol, plo::Array{T}, thres::Array{T}=-plo) where T<:Real
    size(plo) == size(thres) || error("Inconsistent vector lengths")
    tuples = [ber_famiss(r, f, x) for x in plo]
    return [ [ x[y] for x in tuples] for y in 1:2]
end

## actual errors
ber_famiss(r::Roc, plo::ArrayOrReal, thres=-plo) = ber_famiss(r, :θ, plo, thres)

ber(r::Roc, plo::Real, thres::Real=-plo) = +(ber_famiss(r, plo, thres)...)
ber(r::Roc, plo::Real, thres::Array{T}) where {T<:Real} = [ber(r, plo, x) for x in thres]
ber(r::Roc, plo::Array{T}, thres::Array{T}=-plo) where {T<:Real} = [ber(r, x...) for x in zip(plo,thres)]

## minimum ber is best computed using a Roc structure
minber_famiss(r::Roc, plo::ArrayOrReal) = ber_famiss(r, :llr, plo)
"""
`minber()`: The minimum Byes Error Rate for a set of socres, found
after an optimal score-to-llr (log-likelihood-ratio) transformation.
The optimal transformation corresponds to the convex hull of the ROC,
where the optimal llr values correspont to the negative log of the
slope of the appropriate line segment spanning the ROC.  This is
equivalent to running the pool-adjacent-violators algorithm on the
test data.

In order to compute `minber()` for multiple cost functions, as in a Normalized Bayes Error Rate plot,
it is advantageous to compute a `Roc` object once.

Arguments: see `ber()`.
"""
minber(r::Roc, plo::Real) = +(minber_famiss(r, plo)...)
minber(r::Roc, plo::Array{T}) where {T<:Real} = [minber(r, x) for x in plo]

## compute ROC if scores are given
function minber(tar::Vector{T}, non::Vector{T}, plo::ArrayOrReal) where T<:Real
    r = roc(tar, non)
    return minber(r, plo)
end

## factor to normalize the Bayes error rate, for normalized bayes error rate (or normalized cdet)
normfactor(plo) = 1 .+ exp.(abs.(plo))

## factor to convert Bayes error rates to actual / minimum costs.
costfactor(d::DCF) = d.ptar .* d.cmiss + (1 .- d.ptar) .* d.cfa

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
Base.copy(d::DCF) = DCF(d.ptar, d.cfa, d.cmiss)
"""
`setdcf()` allows a global decision cost function to be set.  This DCF can be used in `dcf()` and `mindcf()`.   The decision cost function is defined as

    dcf = ptar Cmiss pmiss + (1-ptar) Cfa pfa.

and is defined by the parameters `ptar` (the target prior), `Cfa` (the cost of a false positive) and `Cmiss` (the cost of a false negative.  Arguments:

 - `ptar`. The target prior, default `0.5`.

 - `cfa`.  The cost of a false positive (false alarm), default `1`

 - `cmiss`.  The cost of a false negative (miss), default `1`.

Multiple simultaneous cost functions can be set by specifying any, or a combination, of these parameters as Vectors.

The current values of the parameters of the DCF can be found using `getdcf()`.
"""
function setdcf(;ptar=0.5, cfa=1, cmiss=1, d=DCF(ptar, cfa, cmiss))
    global global_dcf = copy(d)
end

setdcf(d::DCF) = setdcf(d=d)
## set a default DCF
setdcf()
"""
`getdcf()` retrieves the current valua of the DCF parameters.
"""
getdcf() = global_dcf

"""
`dcf()`.  Compute the classical decision cost function:

    dcf = ptar Cmiss pmiss + (1-ptar) Cfa pfa.

The _relative_ value of the cost function is only determined by the ratio

    oeff = (ptar/(1-ptar) (Cmiss/Cfa).

Therefore, this function is computed using `ber()` and scaling accordingly.  Arguments:

 - `tar`, `non`: target and non-target scores

 - `tnt::TNT`: a `TNT` object containing target and non-target scores.

 - `r::Roc`: a `Roc` object, the result of `roc()`

 - `; d::DCF`: a decision cost function, default `getdcf()`.  This can be a vector or DCFs.

 - `; thres`: the threshold used to make decisions, default `-plo(d)`

 - `; norm`: are the costs normalized to a trivial system deciding using the prior only, default: `false`

"""
dcf(tar::Vector, non::Vector; d::DCF=getdcf(), thres=-plo(d), norm=false) = applyfactor(d, ber(tar, non, plo(d), thres), norm)
dcf(tnt::TNT; d::DCF=getdcf(), thres=-plo(d), norm=false) = dcf(tnt.tar, tnt.non, d=d, thres=thres, norm=norm)
dcf(r::Roc; d::DCF=getdcf(), thres=-plo(d), norm=false) = applyfactor(d, ber(r, plo(d), thres), norm)

"""

`mindcf()` computes the minimum costs of the classical decision cost
function (see `dcf()`) that can be obtained by varying the threshold.
Arguments are the same as for `dcf()`, except that there is no
threshold.

This function uses `minber()` and scales accordingly.
"""
mindcf(r::Roc; d::DCF=getdcf(), norm=false) = applyfactor(d, minber(r, plo(d)), norm)
mindcf(tar::Vector, non::Vector; d::DCF=getdcf(), norm=false) = applyfactor(d, minber(roc(tar, non), plo(d)), norm)
mindcf(tnt::TNT; d::DCF=getdcf(), norm=false) = mindcf(tnt.tar, tnt.non, d=d, norm=norm)

Base.show(io::IO, dcf::DCF) = print(io, "Ptar = ", dcf.ptar, ", Cfa = ", dcf.cfa, ", Cmiss = ", dcf.cmiss)

function Base.show(io::IO, ::MIME"text/plain", dcf::DCF)
    println(io, "Ptar = ", dcf.ptar, ", Cfa = ", dcf.cfa, ", Cmiss = ", dcf.cmiss)
    println(io, " prior log-odds = ", plo(dcf))
    println(io, " effective prior odds = ", oeff(dcf))
    println(io, " effective prior = ", peff(dcf))
end
