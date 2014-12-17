## using Winston

function Winston.plot(r::Roc)
    plot(r.pfa, r.pmiss, "-")
    xlabel("Pfa")
    ylabel("Pmiss")
end

qnorm(x) = √2 * erfinv(2x-1)

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
