## roctypes.jl  Types used for ROC computations
## (c) 2013--2015 David A. van Leeuwen
##
## Licensed under the MIT software license, see LICENSE.md

## TNT: a target-nontarget tuple
"""
    TNT(tar, non)

Holds target and non-target scores, resulting from a test for a two-class classifier.
"""
struct TNT{T<:Real}
    """target scores (scofres for which, actually, the one hypothesis is true)"""
    tar::Vector{T}              # target scores
    """non-target scores (score for which, actually, the oter hypothesis is true)"""
    non::Vector{T}              # non-target scores
end

## Roc: essential data to store Receiver Operating Characteristics
"""
    ROC(pfa::Vector, pmiss::Vector, ch::BitArray, θ::Vector, llr::Vector)

Stores the essential performance information that can be extracted from a
set of supervised trials, i.e., target and non-target scores from a two-class classifier.  Apart from
the (minimalized) arrays for probabilities of false-alarm `pfa` and miss `pmiss`---the coordinates of the ROC
curve---, they are the threshold `θ` for these, a boolean whether-or-not this point lies on the
convex hull `ch`, and the associated optimal log-likelihood-ratio `llr` associated to the line segment.
"""
mutable struct Roc{T<:Real}
    pfa::Vector{Float64}        # probability of false alarm
    pmiss::Vector{Float64}      # probability of miss
    ch::BitArray                # point lies on convex hull
    θ::Vector{T}                # threshold
    llr::Vector{Float64}        # optimal log likelihood ratio
end

## A traditional decision cost function
"""
    DCF(ptar, cfa, cmiss)

Holds one or more cost functions with which the performance of a two-class
classifier can be assessed.  Any of the fields can be either a scalar or an array. The cost function is 
specified by:
  - `ptar`: The prior probability of a target class
  - `cfa`: The cost of a false alarm (false positive)
  - `cmiss`: The cost of a miss (false negative)
"""
struct DCF{PTT,CFT,CMT}
    ptar::PTT               # target prior
    cfa::CFT                # cost of a false alarm
    cmiss::CMT              # cost of a miss
    function DCF{PTT,CFT,CMT}(ptar, cfa, cmiss) where {PTT,CFT,CMT}
        lengths = map(length, (ptar, cfa, cmiss))
        lmax = maximum(lengths)
        for l in lengths
            l in [1,lmax] || error("Inconsistent vector lengths")
        end
        new(ptar, cfa, cmiss)
    end
end

ArrayOrReal{T<:Real} = Union{Array{T}, Real}

DCF(ptar::ArrayOrReal, cfa::ArrayOrReal, cmiss::ArrayOrReal) = DCF{typeof(ptar),typeof(cfa),typeof(cmiss)}(ptar, cfa, cmiss)
