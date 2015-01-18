## ROC.jl  Types and functions for performing ROC analysis
## (c) 2013--2015 David A. van Leeuwen.
##
## Licensed under the MIT software license, see LICENSE.md

module ROC

export ROC, TNT, roc, chllr, eerch, pav, eer, plot, detplot, llrplot, qnorm, pnorm, cllr, mincllr, auc, DataFrame

include("modules.jl")

include("roctypes.jl")

include("includes.jl")

end
