## plot.jl  Placeholders for plotting routines for ROC and DET
## (c) 2014--2015 David A. van Leeuwen
##
## Licensed under the MIT software license, see LICENSE.md

missing(s::AbstractString) = warn(s * "(): you need to load a plotting package (Winston)")

rocplot(r::Roc, nr=1; kwargs...) = missing("rocplot")
detplot(r::Roc, nr=1; kwargs...) = missing("detplot")
llrplot(r::Roc; kwargs...) = missing("llrplot")
apeplot(r::Roc; kwargs...) = missing("apeplot")
nbeplot(r::Roc; kwargs...) = missing("nbeplot")
