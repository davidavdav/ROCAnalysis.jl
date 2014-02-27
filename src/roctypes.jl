## TNT: a target-nontarget tuple
type TNT{T<:Real}
    tar::Vector{T}
    non::Vector{T}
end

type Roc 
    pfa::Vector
    pmiss::Vector
    ch::Vector
    θ::Vector
    llr::Vector
end