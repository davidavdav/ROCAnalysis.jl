## dataframes.jl DataFrames support for ROC
## (c) 2015 David A. van Leeuwen
##
## Licensed under the MIT software license, see LICENSE.md

## We also implement most TNT arguments here.

## Accept a DataFrame with :target and :score columns to most functions
function TNT(x::AbstractDataFrame; score=:score, target=:target)
    t = array(x[target])
    s = array(x[score])
    TNT(s[t], s[!t])
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

dcf(tnt::TNT; d::DCF=getdcf(), thres=-plo(d), norm=false) = dcf(tnt.tar, tnt.non, d=d, thres=thres, norm=norm)
dcf(x::AbstractDataFrame; score=:score, target=:target, d::DCF=getdcf(), thres=-plo(d), norm=false) = dcf(TNT(x, score=score, target=target), d=d, thres=thres, norm=norm)

mindcf(tnt::TNT; d::DCF=getdcf(), norm=false) = mindcf(tnt.tar, tnt.non, d=d, norm=norm)
mindcf(x::AbstractDataFrame; score=:score, target=:target, d::DCF=getdcf(), norm=false) = mindcf(TNT(x, score=score, target=target), d=d, norm=norm)


DataFrames.DataFrame(r::Roc) = DataFrame(pfa=r.pfa, pmiss=r.pmiss, thres=[DataArray(r.θ), NA], chull=r.ch, llr=[DataArray(r.llr), NA])
