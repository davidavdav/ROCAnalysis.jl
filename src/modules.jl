## modules.j. Modules used in ROCAnalysis
## (c) 2015 David A. van Leeuwen.
##
## Licensed under the MIT software license, see LICENSE.md

## Loading modules has been moved to the source files that actually use these modules.

## These are just some functions from NumericalFuns that we somehow had to implement ourselves.

logit(x::AbstractFloat) = log(x / (one(x) - x))
logit(x::Real) = logit(float(x))

softplus(x::AbstractFloat) = x ≤ 0 ? log1p(exp(x)) : x + log1p(exp(-x))
softplus(x::Real) = softplus(float(x))

function sigmoid(x::AbstractFloat)
   if x ≤ 0
	   y = exp(x)
	   return y / (y + one(x))
   else
	   return one(x) / (one(x) + exp(-x))
   end
end
sigmoid(x::Real) = sigmoid(float(x))
