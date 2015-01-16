## dataframes.jl DataFrames support for ROC
## (c) 2015 David A. van Leeuwen
##
## Licensed under the MIT software license, see LICENSE.md

## accept a DataFrame with :target and :score columns to most functions
function TNT(x::DataFrame, score=:score, target=:target)
    t = array(x[target])
    s = array(x[score])
    TNT(s[t], s[!t])
end

for f in (:roc, :eer, :eerch, :auc, :cllr, :mincllr)
    @eval ($f)(tnt::TNT) = ($f)(tnt.tar, tnt.non)
    @eval ($f)(x::DataFrame, score=:score, target=:target) = ($f)(TNT(x, score, target))
end
    
DataFrames.DataFrame(r::Roc) = DataFrame(pfa=r.pfa, pmiss=r.pmiss, thres=[DataArray(r.θ), NA], chull=r.ch, llr=[DataArray(r.llr), NA])
                                       
                                         
