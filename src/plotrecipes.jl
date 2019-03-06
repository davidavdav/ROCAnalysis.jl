using RecipesBase
using Printf

@recipe function plot(r::Roc{T}; convex_hull=false, traditional=false) where T
    ## defaults
    title --> "ROC"
    xlabel --> "Pfa"
    ylabel --> (traditional ? "Phit" : "Pmiss")
    leg --> false
    if convex_hull
        x = r.pfa[r.ch]
        y = r.pmiss[r.ch]
    else
        x = r.pfa
        y = r.pmiss
    end
    if traditional
        y = 1 .- y
    end
    (x, y)
end

@userplot DetPlot

@recipe function detplot(d::DetPlot; convex_hull=false)
    if (length(d.args) != 1 || ! (typeof(d.args[1]) <: Roc))
        error("Argument should be of type Roc, got $(length(d.args)) $(typeof(d.args[1]))", )
    end
    r, = d.args
    ##default
    title --> "DET plot"
    xlabel --> "Pfa (%)"
    ylabel --> "Pmiss (%)"
    leg --> false
    xlim --> (qnorm(0.001), qnorm(0.5))
    ylim --> (qnorm(0.001), qnorm(0.5))
    linewidth --> 2
    ## Not trivial in Plots.jl to control ticks and ticklabels
    ticks --> qnorm.([0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 40] ./ 100)
    formatter --> x -> begin
       p = 100 * pnorm(x)
       if p < 1
           Printf.@sprintf("%3.1f", p)
       else
           Printf.@sprintf("%2.0f", p)
       end
   end
    ## compulsory
    aspect_ratio := :equal

    if convex_hull
        x = qnorm.(r.pfa[r.ch])
        y = qnorm.(r.pmiss[r.ch])
    else
        x = qnorm.(r.pfa)
        y = qnorm.(r.pmiss)
    end
    (x, y)
end

@userplot ApePlot

@recipe function apeplot(ape::ApePlot; xmin = -7, xmax = 7)
    if !(length(ape.args) == 1 && typeof(ape.args[1]) <: Roc)
        error("Expecting one argument of type Roc, got $(typeof(ape.args[1]))")
    end
    r, = ape.args
    title --> "Applied Probability of Error"
    xlabel --> "prior log odds"
    ylabel --> "Bayes' error rate"
    labels --> ["Bayes' Error", "Minimum", "Trivial"]
    linecolors --> [:red :green :black]
    lo = collect(xmin:0.01:xmax)
    be = ber(r, lo)
    mbe = minber(r, lo)
    defbe = 1 ./ ROCAnalysis.normfactor(lo)
    ylim --> (0, 1.1maximum(be))
    x = lo
    y = hcat(be, mbe, defbe)
    (x,y)
end

@userplot NbePlot

@recipe function nbeplot(nbe::NbePlot; xmin = -10, xmax = 5)
    if !(length(nbe.args) == 1 && typeof(nbe.args[1]) <: Roc)
        error("Expecting one argument of type Roc, got $(typeof(nbe.args[1]))")
    end
    r, = nbe.args
    title --> "Normalized Bayes' Error"
    xlabel --> "prior log odds"
    ylabel --> "normalized Bayes' error"
    labels --> ["Bayes' Error", "Minimum BE", "False Alarms", "Misses", "Minimum FA", "Minimum miss"]
    leg --> :left
    linestyles --> [:solid :solid :dash :dot :dash :dot]
    linecolors --> [:red :green :red :red :green :green]
    lo = collect(xmin:0.01:xmax)
    norm = ROCAnalysis.normfactor(lo)
    bfa, bmiss = [be .* norm for be in ROCAnalysis.ber_famiss(r, lo)]
    nbe = bfa + bmiss
    mbfa, mbmiss = [mbe .* norm for mbe in ROCAnalysis.minber_famiss(r, lo)]
    mnbe = minber(r, lo) .* norm
    x = lo
    y = hcat(nbe, mnbe, bfa, bmiss, mbfa, mbmiss)
    (x, y)
end

@userplot LlrPlot

@recipe function llrplot(llr::LlrPlot)
    if !(length(llr.args) == 1 && typeof(llr.args[1]) <: Roc)
        error("Expecting one argument of type Roc, got $(typeof(nbe.args[1]))")
    end
    r, = llr.args
    title --> "LLR plot"
    xlabel --> "score"
    ylabel --> "log likelihood ratio"
    leg --> false
    llr_extrema = extrema(filter(x -> !isinf(x), r.llr))
    mi = max(llr_extrema[1], minimum(r.θ)) ## max of mins
    ma = min(llr_extrema[2], maximum(r.θ)) ## min of maxs
    xlim --> (mi, ma)
    ylim --> (mi, ma)
    (r.θ, r.llr)
end
