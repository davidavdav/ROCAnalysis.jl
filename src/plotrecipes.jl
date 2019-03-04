using RecipesBase
using Printf

@recipe function f(r::Roc{T}; convex_hull=false) where T
    ## defaults
    title --> "ROC"
    xlabel --> "Pfa"
    ylabel --> "Pmiss"
    leg --> false
    if convex_hull
        x := r.pfa[r.ch]
        y := r.pmiss[r.ch]
    else
        x := r.pfa
        y := r.pmiss
    end
    ()
end

## rudimentary DET plot
@userplot DetPlot

@recipe function f(d::DetPlot; convex_hull=false)
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
