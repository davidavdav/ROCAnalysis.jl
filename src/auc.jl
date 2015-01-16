## auc.jl (c) Area Under the ROC
## (c) 2015 David A. van Leeuwen
##
## Licensed under the MIT software license, see LICENSE.md

## The probability that a random non-target score is higher than a target score
auc(r::Roc) = -dot(r.pmiss[1:end-1] + r.pmiss[2:end], diff(r.pfa)) / 2

auc(tar::Vector, non::Vector) = auc(roc(tar, non))
