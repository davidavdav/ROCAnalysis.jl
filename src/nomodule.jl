## nomodul.jl  Separate Module encapsulation from code
## (c) 2013--2015 David A. van Leeuwen
## Licensed under the MIT software license, see LICENSE.md

## re-load this file during development of the package. 

include("modules.jl")

## skip reading types that have been read before
if !isdefined(:TNT)
    include("roctypes.jl")
end

include("includes.jl")
