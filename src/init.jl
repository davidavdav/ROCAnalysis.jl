import Requires

function __init__()::Nothing
    Requires.@require PyPlot="d330b81b-6aea-500a-939a-2ce795aea3ee" include("pyplot.jl")

    Requires.@require Winston="bd07be1c-e76f-5ff0-9c0b-f51ef45303c6" include("winston.jl")

    return nothing
end
