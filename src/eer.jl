## eer.jl (c) 2014 David A. van Leeuwen
## This code is licenced under the MIT software license

## eerch: computes the eer using the convex hull method
eerch(r::Roc) = eerch(r.pfa, r.pmiss, r.ch)

## compute the EER given corresponding pfa and pmiss array, but use the
## cunvex hull points only.
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

## just a simple EER approximation, optimized for memory and speed
function eer{T}(tar::Vector{T}, non::Vector{T})
    ntar = length(tar)
    nnon = length(non)
    scores = vcat(tar,non)      # targets before non-targets
    so = sortperm(scores)
    eer(so, ntar, ntar)
end

function eer{T<:Integer}(so::Vector{T}, ntar::T, selection::BitVector)
    length(so) == length(selection) || error("selection not the same length as sort order")
    eer(so[findin(so, find(selection))], ntar)
end

## This computes the actual approximate eer from sortorder and number of targets and nontargets
## The argument "so" may be a subset of an original sort order, "thres" should be the
## original number of targets
function eer{T<:Integer}(so::Vector{T}, thres::T, ntar=sum(so .<= thres))
    nnon = length(so) - ntar
    Δfa = 1. / nnon
    Δmiss = 1. / ntar
    pfa = minerr = mindiff = 1.
    pmiss = absΔ = 0.
    for s in so
        if s <= thres            # target score
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

## Returns the index to the largest value in the sorted array `a` ≤ `x` if !lower
## If lower==true, the value must be strictly < `x`
function binsearch{T}(x::T, a::Vector{T}, lower=false)
    issorted(a) || error("Array needs to be sorted")
    mi = 1
    ma = length(a)
    if x < a[mi] || lower && x == a[mi]
        return 0
    elseif x > a[ma] || !lower && x == a[ma]
        return ma
    end
    while ma - mi > 1
        h = mi + div(ma-mi,2)
        if x > a[h] || !lower && x == a[h]
            mi = h
        else
            ma = h
        end
    end
    return mi
end

## if tar and non are sorted, we may be doing even faster...
function eer_sorted{T}(tar::Vector{T}, non::Vector{T})
    ntar = length(tar)
    nnon = length(non)
    Δfa = 1. / nnon
    Δmiss = 1. / ntar
    absΔ = mindiff = minerr = 1.
    pmiss = 0.
    ni = binsearch(tar[1], non, true)
    pfa = 1. - ni * Δfa
    ti = 1
    ## non[ni] < tar[ti] ≤ non[ni+1]
    ni += 1
    while ni ≤ nnon && ti ≤ ntar
        if tar[ti] < non[ni]
            pmiss += Δmiss
            ti += 1
        elseif tar[ti] > non[ni]
            pfa -= Δfa
            ni += 1
        else
            pmiss += Δmiss
            ti += 1
            pfa -= Δfa
            ni += 1
        end
        absΔ = abs(pfa - pmiss)
        if absΔ < mindiff
            mindiff = absΔ
            minerr = (pfa + pmiss) / 2
        end
    end
    minerr
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

