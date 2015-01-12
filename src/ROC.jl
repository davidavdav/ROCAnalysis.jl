## ROC.jl  Types and functions for performing ROC analysis
## (c) 2013 David A. van Leeuwen. 

module ROC

export ROC, TNT, roc, chllr, eerch, pav, eer, detplot, llrplot, qnorm, pnorm, cllr, minclr, auc

using CHull
using Winston
using NumericExtensions

include("roctypes.jl")
include("roc.jl")
include("eer.jl")
include("cllr.jl")
include("auc.jl")
include("pav.jl")
include("plot.jl")

end
