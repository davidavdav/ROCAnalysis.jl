## roc.jl  Receiver Operating Characteristics computations
## (c) 2013--2015 David A. van Leeuwen
##
## Licensed under the MIT software license, see LICENSE.md
## Routines to compute Receiver Operating Characteristics

## These routines need a re-write, we can probably unroll a few loops
## and be more efficient with memory and sorting.

function roc{T<:Real}(tar::Vector{T}, non::Vector{T}; laplace::Bool=false) 
    xo, tc = sortscores(tar, non)
    ## first collect point of the same threshold (for discrete score data)
    θ, tc, nc, = binscores(xo, tc)
    ## now we can compute pmiss and pfa
    pfa = [1, 1 - cumsum(nc)/length(non)]
    pmiss = [0, cumsum(tc)/length(tar)]
    ## convex hull and llr
    ch, llr = chllr(tc, nc, xo, laplace=laplace)
    ## remove points on the ROC in the middle of a horizontal or vertical segment
    changes = changepoints(pfa, pmiss)
    (pfa, pmiss, ch) = map(a -> a[changes], (pfa, pmiss, ch))
    (llr, θ) = map(a -> a[changes[1:end-1]], (llr, θ))
    Roc(pfa, pmiss, ch, θ, llr)
end

## this returns the target count (1's an 0's for targets and non targets) and scores, 
## in the order of the scores. 
function sortscores{T<:Real}(tar::Vector{T}, non::Vector{T})
    x = vcat(tar, non)                     # order targets before non-targets
    so = sortperm(x)                       # sort order
    xo = x[so]                             # scores ordered
    tc = vcat(ones(Int, length(tar)), zeros(Int, length(non)))[so] # target count
    return xo, tc
end    

## bins scores that are the same into aggredated score counts
function binscores{T<:Real}(xo::Vector{T},tc::Vector{Int})
    nc = 1 - tc
    changes = find([true, diff(xo) .!= 0, true]) # points where threshold change
    θ = xo                                       # threshold
    keep = trues(length(tc))
    @inbounds for i=1:length(changes)-1
        start = changes[i]
        stop = changes[i+1]-1
#        for j=start+1:stop
#            tc[start] += tc[j]
#            nc[start] += nc[j]
#            keep[j] = false
#        end
        if stop>start
            tc[start] = sum(tc[start:stop])
            nc[start] = sum(nc[start:stop])
            keep[start+1:stop]= false
        end
    end
    (tc, nc, θ) = map(a -> a[keep], (tc, nc, θ))
    return θ, tc, nc, keep
end
## This finds the change points on the ROC (the corner points)
function changepoints(pfa::Vector{Float64}, pmiss::Vector{Float64})
    (const_pfa, const_pmiss) = map(a -> diff(a) .== 0, (pfa, pmiss)) # no changes
    return [true, ! (const_pfa[1:end-1]  &  const_pfa[2:end] | 
                     const_pmiss[1:end-1] & const_pmiss[2:end]), true]
end

## compute convex hull and optimal llr from target count, nontarget count, and ordered scores
function chllr{T}(tc::Vector{Int}, nc::Vector{Int}, xo::Vector{T}; laplace::Bool=true)
    if laplace
        tc = [1, 0, tc, 1, 0]
        nc = [0, 1, nc, 0, 1]
        xo = [-Inf, -Inf, xo, Inf, Inf]
    end
    ntar = sum(tc)
    nnon = sum(nc)
#    if fast
        pfa = [1, 1 - cumsum(nc)/nnon] # use rationals for accuracy of the convex hull
        pmiss = [0, cumsum(tc)/ntar]
        index = rochull(pfa, pmiss)
#    else
#        pfa = [1, 1 - cumsum(nc)/nnon]
#        pmiss = [0, cumsum(tc)/ntar]
        ## convex hull
#        hull = chull(vcat(hcat(pfa, pmiss), [2 2])) # convex hull points
#        index = sort(hull.vertices[:,1])[1:end-1]   # indices of the points on the CH
#    end
    ch = falses(length(pfa))
    ch[index] = true
    ## LLR
    (Δpfa, Δpmiss) = map(a -> diff(a[ch]), (pfa, pmiss))
    llr = log(abs(Δpmiss ./ Δpfa))[cumsum(ch)[1:end-1]]
    if laplace
        return [true, ch[4:end-3], true], llr[3:end-2]
    else
        return ch, llr
    end
end

## returns positive number if test point is "left" of line (x1,y1) -- (x2,y2)
## positive: left, negative: right, zero: on
## This works best with rationals, but that makes it slow 
function isleft{T<:Real}(x1::T, y1::T, x2::T, y2::T, xt::T, yt::T)
    (x2 - x1)*(yt - y1) - (xt - x1)*(y2 - y1) 
end

## Andrew's Monotone Chain Algorithm, from http://geomalgorithms.com/a10-_hull-1.html
## (ignoring the C-code:-)
## points are (pmiss, pfa), pmiss (x) is increasing, pfa (y) is decreasing
## we implictly add a point at (2,2) that is on the convex hull
function rochull{T<:Real}(pfa::Vector{T}, pmiss::Vector{T})
    two = 2one(T)
    n = length(pmiss)
    ## minmin and minmax
    i = 1
    while i < n
        if pmiss[i+1] == pmiss[1]
            i += 1
        else
            break
        end
    end
    minmin = i
    minmax = 1
    ## maxmin and maxmax
    maxmin = maxmax = n+1       # virtual point for now
    stack = [minmin]
    for i=minmin+1:n
        if isleft(pmiss[minmin], pfa[minmin], two, two, pmiss[i], pfa[i]) ≥ zero(T)
            continue
        end
        while length(stack) ≥ 2
            if isleft(pmiss[stack[end-1]], pfa[stack[end-1]],
                      pmiss[stack[end]], pfa[stack[end]], pmiss[i], pfa[i]) > zero(T)
                break
            else
                pop!(stack)
            end
        end
        push!(stack, i)
    end
    unshift!(stack, minmax)
    stack
end

