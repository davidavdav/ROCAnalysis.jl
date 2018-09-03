## winston.jl  Winston plotting routines for ROC and DET
## (c) 2014--2015 David A. van Leeuwen
##
## Licensed under the MIT software license, see LICENSE.md

import Winston

function rocplot(r::Roc, nr=1; ch=false, title="ROC", xlabel="Pfa", ylabel="Pmiss")
    col = string("krgmcb"[(nr-1) % 6 + 1])
    lty = ["-", "--", ";"][div(nr-1,6) + 1]
    if ch
        chi = find(r.ch)
        x, y = r.pfa[chi], r.pmiss[chi]
    else
        x, y = r.pfa, r.pmiss
    end
    if nr==1
        p = Winston.plot(x, y, lty*col)
        Winston.xlabel(xlabel)
        Winston.ylabel(ylabel)
        Winston.title(title)
    else
        p = Winston.oplot(x, y, lty * col)
    end
    return p
end
Winston.plot(r::Roc, nr=1; ch=false, title="ROC", xlabel="Pfa", ylabel="Pmiss") = rocplot(r::Roc, nr; title=title, xlabel=xlabel, ylabel=ylabel)

function detplot(r::Roc, nr=1; xlabel="Pfa (%)", ylabel="Pmiss (%)", title="DET plot")
    if nr==1
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
        p = Winston.plot()
        Winston.setattr(p, :aspect_ratio, 1) # essential for DET plots!
        Winston.title(title)
        Winston.xlim(qnorm(0.8e-3), qnorm(0.55))
        Winston.ylim(qnorm(0.8e-3), qnorm(0.55))
        for axis in (p.x1, p.y1)
            Winston.setattr(axis, "draw_subticks", false)
            Winston.setattr(axis, "ticks", qnorm(grid))
            Winston.setattr(axis, "ticklabels", ticklabels)
            Winston.setattr(axis, "tickdir", 1)
            Winston.setattr(axis, "draw_grid", true)
        end
        for axis in (p.x2, p.y2)
            Winston.setattr(axis, "draw_ticks", false)
        end
        Winston.xlabel(xlabel)
        Winston.ylabel(ylabel)
    end
    col = string("krgmcb"[(nr-1) % 6 + 1])
    lty = ["-", "--", ";"][div(nr-1,6) + 1]
    p = Winston.oplot(qnorm(r.pfa), qnorm(r.pmiss), lty * col)
    return p
end

function llrplot(r::Roc; title="LLR plot", xlabel="score", ylabel="LLR")
    mi = max(minimum(r.llr), minimum(r.θ))
    ma = min(maximum(r.llr), maximum(r.θ))
    ran = [mi, ma]
    p = Winston.plot(r.θ, r.llr)
    Winston.xlim(ran)
    Winston.ylim(ran)
    Winston.title(title)
    Winston.xlabel(xlabel)
    Winston.ylabel(ylabel)
    return p
end

function apeplot(r::Roc; xmin=-7, xmax=7, title="Applied Probability of Error", xlabel="prior log odds", ylabel="Bayes error rate")
    lo = collect(xmin:0.01:xmax)
    be = ber(r, lo)
    mbe = minber(r, lo)
    defbe = 1 ./ normfactor(lo)
    p = Winston.plot(lo, be, "-r")
    Winston.oplot(lo, mbe, "-g")
    Winston.oplot(lo, defbe, "--k")
    Winston.ylim(0, 1.1*maximum(be))
    Winston.title(title)
    Winston.xlabel(xlabel)
    Winston.ylabel(ylabel)
    return p
end

function nbeplot(r::Roc; xmin = -10, xmax = 5, title="Normalized Bayes' Error", xlabel="prior log odds", ylabel="normalized Bayes error")
    lo = collect(xmin:0.01:xmax)
    norm = normfactor(lo)
    bfa, bmiss = [be .* norm for be in ber_famiss(r, lo)]
    nbe = bfa + bmiss
    mbfa, mbmiss = [mbe .* norm for mbe in minber_famiss(r, lo)]
    mnbe = minber(r, lo) .* norm
    p = Winston.plot(lo, nbe, "-r")
    Winston.oplot(lo, mnbe, "-g")
    Winston.oplot(lo, bfa, "--r")
    Winston.oplot(lo, bmiss, "-.r")
    Winston.oplot(lo, mbfa, "--g")
    Winston.oplot(lo, mbmiss, "-.g")
    Winston.ylim(0, 1.1)
    Winston.title(title)
    Winston.xlabel(xlabel)
    Winston.ylabel(ylabel)
    return p
end
