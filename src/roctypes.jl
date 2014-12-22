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
