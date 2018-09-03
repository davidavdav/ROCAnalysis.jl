## eer.jl  Equal Error Rate and friends
## (c) 2014--2015 David A. van Leeuwen
##
## This code is licenced under the MIT software license

## eerch: computes the eer using the convex hull method
"""
`eerch()` computes the Equal Error Rate using the Convex Hull interpretation.  This computes the error
rate at which the convex hull of the ROC curve crosses the `y=x` diagonal in the ROC.

This value is less sensitive to trial sampling at low number of trials, and has the interpretation of
`the maximum while vary prior of the minimum decision costs'---a useful quantity in decision cost
function analysis of a two class classification system.

Arguments:

 - `r::Roc`: an object of type `Roc`, the result of `roc()`

 - `tar, non`: Vectors of target and non-target scores

 - `pfa, pmiss, ch`: Vectors of the false positive and false negative rate, and an array indicating
  the membership of the `(pfa, pmiss)` point on the convex hull.  These points are found by `roc()`.
"""
eerch(r::Roc) = eerch(r.pfa, r.pmiss, r.ch)
eerch(tar::Vector, non::Vector) = eerch(roc(tar, non))

## compute the EER given corresponding pfa and pmiss array, but use the
## convex hull points only.
function eerch(pfa::Vector{T}, pmiss::Vector{T}, ch::BitVector) where T<:AbstractFloat
    @assert length(pfa) == length(pmiss) == length(ch)
    ## find the index on the ch where the crossing of pfa=pmiss occurs
    chi = findall(ch)
    i = chi[1]
    direction = sign(pfa[i] - pmiss[i]) # >0 if first pfa > pmiss
    ## does it pay off to do a binary search?  It should, really.
    li = i
    for loop_variable in 1:length(chi)
        i = chi[loop_variable]
        if sign(pfa[i] - pmiss[i]) != direction
            break
        end
        li = i
    end
    ## compute the crossing
    (ax,bx) = pfa[[li,i]]
    (ay,by) = pmiss[[li,i]]
    return crossing(ax, ay, bx, by)
end

## compute crossing with y=x points (ax,ay) and (bx,by)
function crossing(ax, ay, bx, by)
    den = ax-ay-bx+by
    if den == 0
        return (ax + bx + ay + by) / 4
    else
        return ax + (ax-ay)*(bx-ax) / den
    end
end

## just a simple EER approximation, optimized for memory and speed
"""
`eer(tar, non)` computes a simple approximation to the Equal Error
Rate, the intersection of the ROC with the line `y=x`, i.e., the error
rate at which the false positive rate and the false negative rate are
approximately equal.

For a more consistent interpretation and implementation of the EER, please consider `eerch()`,
the Equal Error Rate - Convex Hull computation.
"""
eer(tar::Vector{T}, non::Vector{T}) where {T<:Real} = eer_sorted(sort(tar), sort(non))
eer(r::Roc) = eerch(r.pfa, r.pmiss, trues(length(r.pfa)))

## if tar and non are sorted, we may be doing even faster...
function eer_sorted(tar::Vector{T}, non::Vector{T}) where T<:Real
    ntar = length(tar)
    nnon = length(non)
    ntar > 0 && nnon > 0 || error("target and non-target arrays must not be empty")
    Δfa = 1. / nnon
    Δmiss = 1. / ntar
    pmiss = 0.
    ni = binsearch(tar[1], non, lower=true) # index in non-targets of lowest target score
    pfa = 1. - ni * Δfa
    ti = 1
    ni += 1
    lpfa, lpmiss = pfa, pmiss
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
        if pfa < pmiss
            break
        end
        lpfa, lpmiss = pfa, pmiss
    end
    return crossing(pfa, pmiss, lpfa, lpmiss)
end

## Returns the index to the largest value in the sorted array `a` ≤ `x` if lower==false
## If lower==true, the value must be strictly < `x`
function binsearch(x::Real, a::Vector{T}; lower=false, check=false) where T<:Real
    check && (issorted(a) || error("Array needs to be sorted"))
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

## EER, given the sort order of `[tar, non]`
## The idea is that you don't need to sort multiple times if you have
## multiple selections.  This cannot be accurate for cases where
## `tar` and `non` have the same value
function eer_so(so::Vector{T}, ntar::T, selection::BitVector) where T<:Integer
    length(so) == length(selection) || error("selection not the same length as sort order")
    eer(so[findin(so, find(selection))], ntar)
end

## This computes the actual approximate eer from sortorder and number of targets and nontargets
## The sort order "so" may be a subset of an original sort order, "thres" should be the
## original number of targets.
function eer_so(so::Vector{T}, thres::T, ntar::T=sum(so .<= thres)) where T<:Integer
    nnon = length(so) - ntar
    Δfa = 1. / nnon
    Δmiss = 1. / ntar
    pfa = 1.
    pmiss = 0.
    lpfa, lpmiss = pfa, pmiss
    for s in so
        if s <= thres            # target score
            pmiss += Δmiss
        else
            pfa -= Δfa
        end
        if pfa < pmiss
            break
        end
        lpfa, lpmiss = pfa, pmiss
    end
    return crossing(pfa, pmiss, lpfa, lpmiss)
end

## experimental routine for timing and optimization...
function eer_experimental(tar::Vector{T}, non::Vector{T}; method=:naive, fast=true) where T
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
