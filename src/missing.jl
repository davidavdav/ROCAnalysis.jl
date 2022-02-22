## missing.jl (c) Area Under the ROC
## (c) 2016 David A. van Leeuwen
##
## Licensed under the MIT software license, see LICENSE.md

remove_missing(array::AbstractVector{Union{T,Missing}}) where T = T[value for value in array if ! isa(value, Missing)]
