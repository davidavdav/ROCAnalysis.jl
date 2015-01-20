## plot.jl  Plotting routines for ROC and DET
## (c) 2014--2015 David A. van Leeuwen
##
## Licensed under the MIT software license, see LICENSE.md

function Winston.plot(r::Roc, nr=1; ch=false)
    col = string("krgmcb"[(nr-1) % 6 + 1])
    lty = ["-", "--", ";"][div(nr-1,6) + 1]
    if ch
        chi = find(r.ch)
        x, y = r.pfa[chi], r.pmiss[chi]
    else
        x, y = r.pfa, r.pmiss
    end
    if nr==1
        plot(x, y, lty*col)
        xlabel("Pfa")
        ylabel("Pmiss")
    else
        oplot(x, y, lty * col)
    end
end

## R terminology, quantile function qnorm() and cumulative distribution pnorm()
qnorm(x) = √2 * erfinv(2x-1)
pnorm(x) = (1 + erf(x/√2)) / 2 

function detplot(r::Roc, nr=1)
    if nr==1
        grid = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 40]
        ticklabels = String[]
        for x in grid
            if x<1
                push!(ticklabels, @sprintf("%3.1f", x))
            else
                push!(ticklabels, @sprintf("%1.0f", x))
            end
        end
        grid ./= 100
        p = plot()
        title("DET plot")
        xlim(qnorm(0.8e-3), qnorm(0.55))
        ylim(qnorm(0.8e-3), qnorm(0.55))
        for axis in (p.x1, p.y1)
            setattr(axis, "draw_subticks", false)
            setattr(axis, "ticks", qnorm(grid))
            setattr(axis, "ticklabels", ticklabels)
            setattr(axis, "tickdir", 1)
            setattr(axis, "draw_grid", true)
        end
        for axis in (p.x2, p.y2)
            setattr(axis, "draw_ticks", false)
        end
        xlabel("Pfa (%)")
        ylabel("Pmiss (%)")
    end
    col = string("krgmcb"[(nr-1) % 6 + 1])
    lty = ["-", "--", ";"][div(nr-1,6) + 1]
    oplot(qnorm(r.pfa), qnorm(r.pmiss), lty * col)
end

function llrplot(r::Roc)
    mi = max(minimum(r.llr), minimum(r.θ))
    ma = min(maximum(r.llr), maximum(r.θ))
    ran = [mi,ma]
    plot(r.θ, r.llr)
    xlim(ran)
    ylim(ran)
    title("LLR plot")
    xlabel("score")
    ylabel("LLR")
end

function apeplot(r::Roc; xmin=-7, xmax=7)
    lo = [xmin:0.01:xmax]
    be = ber(r, lo)
    mbe = minber(r, lo)
    defbe = 1 ./ normfactor(lo)
    plot(lo, be, "-r")
    oplot(lo, mbe, "-g")
    oplot(lo, defbe, "--k")
    ylim(0, 1.1*maximum(be))
    title("Applied Probability of Error")
    xlabel("prior log odds")
    ylabel("Bayes error rate")
end

function nbeplot(r::Roc; xmin = -10, xmax = 5)
    lo = [xmin:0.01:xmax]
    norm = normfactor(lo)
    bfa, bmiss = [be .* norm for be in ber_famiss(r, lo)]
    nbe = bfa + bmiss
    mbfa, mbmiss = [mbe .* norm for mbe in minber_famiss(r, lo)]
    mnbe = minber(r, lo) .* norm
    plot(lo, nbe, "-r")
    oplot(lo, mnbe, "-g")
    oplot(lo, bfa, "--r")
    oplot(lo, bmiss, "-.r")
    oplot(lo, mbfa, "--g")
    oplot(lo, mbmiss, "-.g")
    ylim(0, 1.1)
    title("Normalized Bayes' Error")
    xlabel("prior log odds")
    ylabel("normalized Bayes error")
end
