## ROCAnalysis.jl  Types and functions for performing ROC analysis
## (c) 2013--2016 David A. van Leeuwen.
##
## Licensed under the MIT software license, see LICENSE.md

module ROCAnalysis

export Roc, DET, TNT, DCF, DataFrame, roc, chllr, eerch, pav, eer, plot, detplot, llrplot, apeplot, nbeplot, qnorm, pnorm, cllr, mincllr, auc, AUC, oeff, peff, plo, ber, minber, dcf, mindcf, setdcf, getdcf, sigmoid

include("modules.jl")

include("roctypes.jl")

include("includes.jl")

end
