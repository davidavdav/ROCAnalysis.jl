## pyplot.jl  PyPlot plotting routines for ROC and DET
## (c) 2016 David A. van Leeuwen
##
## Licensed under the MIT software license, see LICENSE.md

import PyPlot

function detplot(r::Roc, nr=1; xlabel="Pfa (%)", ylabel="Pmiss (%)", title="DET plot")
    if nr==1
        fig = PyPlot.figure("DET plot")
        PyPlot.title("DET plot (must fix axes)")
        grid = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 40]
        ticklabels = AbstractString[]
        for x in grid
            if x<1
                push!(ticklabels, @sprintf("%3.1f", x))
            else
                push!(ticklabels, @sprintf("%1.0f", x))
            end
        end
        grid ./= 100
        PyPlot.axes(aspect="equal")
        ## This should be done with matplotlib scale, but I don't know the Julia interface.
        PyPlot.axis([qnorm(0.8e-3), qnorm(0.55), qnorm(0.8e-3), qnorm(0.55)])
        PyPlot.plot(qnorm(r.pfa), qnorm(r.pmiss))
        PyPlot.xlabel(xlabel)
        PyPlot.ylabel(ylabel)

    end
end
