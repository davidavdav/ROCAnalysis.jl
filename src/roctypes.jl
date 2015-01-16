## roctypes.jl  Types used for ROC computations
## (c) 2013--2015 David A. van Leeuwen
##
## Licensed under the MIT software license, see LICENSE.md


## TNT: a target-nontarget tuple
type TNT{T<:Real}
    tar::Vector{T}
    non::Vector{T}
end

type Roc{T<:Real}
    pfa::Vector{Float64}
    pmiss::Vector{Float64}
    ch::BitArray
    θ::Vector{T}
    llr::Vector{Float64}
end
