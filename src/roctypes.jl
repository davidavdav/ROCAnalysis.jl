## TNT: a target-nontarget tuple
type TNT{T<:Real}
    tar::Vector{T}
    non::Vector{T}
end

type Roc{T<:Real}
    pfa::Vector{T}
    pmiss::Vector{T}
    ch::BitArray
    θ::Vector{T}
    llr::Vector{T}
end
