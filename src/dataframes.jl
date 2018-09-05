## dataframes.jl DataFrames support for ROC
## (c) 2015 David A. van Leeuwen
##
## Licensed under the MIT software license, see LICENSE.md

using DataFrames: AbstractDataFrame, DataFrame

## We also implement most TNT arguments here.

## Accept a DataFrame with :target and :score columns to most functions
"""
TNT(::DataFrame; optional-args) extracts target and non-target scores from a dataframe.  The optional arhuments encode which colums are used for what purpose. `score` (default `:score`) is the name of the column containing the classifier's scores.  `target` (default `:target`), of type `Bool`, determines whether this is a target score or not.
"""
function TNT(x::AbstractDataFrame; score=:score, target=:target)
    ok = .!(ismissing.(x[target]) .| ismissing.(x[score]))
    t = BitArray(x[target][ok])
    s = collect(skipmissing(x[score][ok]))
    TNT(s[t], s[.!t])
end

for f in (:eer, :eerch, :auc, :cllr, :mincllr)
    @eval ($f)(tnt::TNT) = ($f)(tnt.tar, tnt.non)
    @eval ($f)(x::AbstractDataFrame; score=:score, target=:target) = ($f)(TNT(x, score=score, target=target))
end

roc(x::AbstractDataFrame; score=:score, target=:target, laplace=false, collapse=true) = roc(TNT(x, score=score, target=target), laplace=laplace, collapse=collapse)
roc(tnt::TNT; laplace=false, collapse=true) = roc(tnt.tar, tnt.non, laplace=laplace, collapse=collapse)

for f in (:ber, :minber)
    @eval ($f)(tnt::TNT, plo) = ($f)(tnt.tar, tnt.non, plo)
    @eval ($f)(x::AbstractDataFrame, plo; score=:score, target=:target) = ($f)(TNT(x, score=score, target=target), plo)
end

dcf(x::AbstractDataFrame; score=:score, target=:target, d::DCF=getdcf(), thres=-plo(d), norm=false) = dcf(TNT(x, score=score, target=target), d=d, thres=thres, norm=norm)

mindcf(x::AbstractDataFrame; score=:score, target=:target, d::DCF=getdcf(), norm=false) = mindcf(TNT(x, score=score, target=target), d=d, norm=norm)

import DataFrames.DataFrame
"""
`Dataframe(::Roc)` converts a `Roc` object in a datframe with columns:

 - `pfa` the false alarm rate

 - `pmiss` the miss rate

 - `thres` the thereshold, separating this line's pfa and pmiss from the next

 - `chull` indicating if this point is on the convex hull of the ROC curve

 - `llr` the optimal log-likelihood-ratio score for all data points contributing to the ROC line segment from this line to the next
"""
DataFrame(r::Roc) = DataFrame(pfa=r.pfa, pmiss=r.pmiss, thres=[Array(r.θ); missing], chull=r.ch, llr=[Array(r.llr); missing])
