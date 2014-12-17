## eerch: computes the eer using the convex hull method
function eerch{T<:FloatingPoint}(pfa::Vector{T}, pmiss::Vector{T}, ch::BitVector)
    @assert length(pfa) == length(pmiss) == length(ch)
    ## find the index on the ch where the crossing of pfa=pmiss occurs
    i = find(diff(sign(pfa[ch] - pmiss[ch])) .!= 0)[1]
    if i==length(ch)
        warning("EER at the extreme")
        return (pfa[i] + pmiss[i])/2
    end
    ## compute the crossing
    (ax,bx) = pfa[ch][i:i+1]
    (ay,by) = pmiss[ch][i:i+1]
    return ax + (ax-ay)*(bx-ax) / (ax-ay-bx+by)
end

## just a simple EER approximation, optimized for memory and speech
function eer{T}(tar::Vector{T}, non::Vector{T})
    ntar = length(tar)
    nnon = length(non)
    scores = vcat(tar,non)      # targets before non-targets
    so = sortperm(scores)
    eer(T, so, ntar, ntar)
end

## This computes the actual approximate eer from sortorder and number of targets and nontargets
## The argument "so" may be a subset of an original sort order, "thres" should be the
## original number of targets
function eer(T::Type, so::Vector, thres::Int, ntar=sum(so .<= thres))
    nnon = length(so) - ntar
    Δfa = one(T) / nnon
    Δmiss = one(T) / ntar
    minerr = mindiff = one(T)
    pfa = one(T)
    pmiss = zero(T)
    absΔ = zero(T)
    for s in so
        if s <= ntar            # target score
            pmiss += Δmiss
        else
            pfa -= Δfa
        end        
        absΔ = abs(pfa - pmiss)
        if absΔ < mindiff
            mindiff = absΔ
            minerr = (pfa + pmiss) / 2
        end
    end
    return minerr
end

function eerstats{T}(tar::Vector{T}, non::Vector{T}; method=:naive, fast=true)
    ntar = length(tar)
    nnon = length(non)
    @time scores = vcat(non,tar)      # 1
    @time so = sortperm(scores)       # 1
    @time truth = so .> nnon          # 1/64
    @time if !fast
        pfa = 1.0 - cumsum(!truth)/nnon
        pmiss = cumsum(truth)/ntar
    else                        # this takes less memory, but that's all I can say for it. 
        Δ = 1/nnon
        pfa = similar(scores)   # 1
        s = one(T)
        for i=1:length(scores)
            s -= Δ*!truth[i]
            pfa[i] = s
        end
        Δ = 1/ntar
        pmiss = similar(scores) # 1
        s = zero(T)
        for i=1:length(scores)
            s += Δ*truth[i]
            pmiss[i] = s
        end
    end
    @time if method==:naive
        i = indmin(abs(pfa-pmiss))
        res = (pfa[i] + pmiss[i]) / 2
    end
    return res, (pfa, pmiss, scores[so]) # 1
end

