## dataframes.jl DataFrames support for ROC
## (c) 2015 David A. van Leeuwen
##
## Licensed under the MIT software license, see LICENSE.md

## accept a DataFrame with :target and :score columns to most functions
function TNT(x::DataFrame; score=:score, target=:target)
    t = array(x[target])
    s = array(x[score])
    TNT(s[t], s[!t])
end

for f in (:eer, :eerch, :auc, :cllr, :mincllr)
    @eval ($f)(tnt::TNT) = ($f)(tnt.tar, tnt.non)
    @eval ($f)(x::DataFrame; score=:score, target=:target) = ($f)(TNT(x, score=score, target=target))
end

roc(x::DataFrame; score=:score, target=:target, laplace=false, collapse=false) = roc(TNT(x, score=score, target=target), laplace=laplace, collapse=collapse)
roc(tnt::TNT; laplace=false, collapse=false) = roc(tnt.tar, tnt.non, laplace=laplace, collapse=collapse)

for f in (:ber, :minber)
    @eval ($f)(tnt::TNT, plo) = ($f)(tnt.tar, tnt.non, plo)
    @eval ($f)(x::DataFrame, plo; score=:score, target=:target) = ($f)(TNT(x, score=score, target=target), plo)
end

for f in (:dcf, :mindcf)
    @eval ($f)(tnt::TNT, d::DCF; norm=false) = ($f)(tnt.tar, tnt.non, d, norm=norm)
    @eval ($f)(x::DataFrame, d::DCF; score=:score, target=:target, norm=false) = ($f)(TNT(x, score=score, target=target), d, norm=norm)
end

DataFrames.DataFrame(r::Roc) = DataFrame(pfa=r.pfa, pmiss=r.pmiss, thres=[DataArray(r.θ), NA], chull=r.ch, llr=[DataArray(r.llr), NA])
