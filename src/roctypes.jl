## roctypes.jl  Types used for ROC computations
## (c) 2013--2015 David A. van Leeuwen
##
## Licensed under the MIT software license, see LICENSE.md


## TNT: a target-nontarget tuple
type TNT{T<:Real}
    tar::Vector{T}              # target scores
    non::Vector{T}              # non-target scores
end

## Roc: essential data to store Receiver Operating Characteristics
type Roc{T<:Real}
    pfa::Vector{Float64}        # probability of false alarm
    pmiss::Vector{Float64}      # probability of miss
    ch::BitArray                # point lies on convex hull
    θ::Vector{T}                # threshold
    llr::Vector{Float64}        # optimal log likelihood ratio
end

## A traditional decision cost function
type DCF{PTT,CFT,CMT}
    ptar::PTT               # target prior
    cfa::CFT                # cost of a false alarm
    cmiss::CMT              # cost of a miss
    function DCF(ptar, cfa, cmiss)
        lengths = map(length, (ptar, cfa, cmiss))
        lmax = maximum(lengths)
        for l in lengths
            l in [1,lmax] || error("Inconsistent vector lengths")
        end
        new(ptar, cfa, cmiss)
    end
end
typealias ArrayOrReal{T<:Real} Union(Array{T}, Real)
DCF(ptar::ArrayOrReal, cfa::ArrayOrReal, cmiss::ArrayOrReal) = DCF{typeof(ptar),typeof(cfa),typeof(cmiss)}(ptar, cfa, cmiss)
