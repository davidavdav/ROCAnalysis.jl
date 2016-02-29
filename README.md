ROCAnalysis.jl
======
[![Build Status](https://travis-ci.org/davidavdav/ROCAnalysis.jl.svg?branch=master)](https://travis-ci.org/davidavdav/ROCAnalysis.jl)

Install
-------

```julia
Pkg.add("ROCAnalysis")
```

Introduction
------------

Receiver Operating Characteristic Analysis functions for evaluation probabilistic binary classifiers.   Allows efficient plotting of ROC curves, and many more things. 

Please note there is an alternative implementation under [a similar name](https://github.com/diegozea/ROC.jl), and support for ROC analysis also exists in [MLBase](https://github.com/lindahua/MLBase.jl). 

Our implementation is more geared towards:
 - large amounts of data, with efficient ROC statistics calculation
 - Decision Cost Function analysis
 - Detection Error Trade-off (DET) analysis
 - ROC convex hull computation, analysis and EER-interpretation
 - Optimal Likelihood Ratio computation

The development roadmap is largely based on the functionality in a similar R package [ROC](https://github.com/davidavdav/ROC). 

Synopsis
----
An [annotated jupyter notebook](ROCAnalysis.ipynb) of the code below can be found [here](http://nbviewer.jupyter.org/github/davidavdav/ROCAnalysis.jl/blob/master/ROCAnalysis.ipynb).
```julia
using ROCAnalysis
## Produce some well-calibrated log-likelihood-ratio scores for target and non-target class:
tar =  2 + 2randn(1000)
non = -2 + 2randn(100000)
## quick estimate of equal error rate, should be close to pnorm(-1) = 0.5 + 0.5erf(-1/√2)
eer(tar, non) 
## compute full ROC statistics
r = roc(tar, non)
## accurate computation of the equal error rate, using the convex hull
eerch(r)
## roc plot, we plot errors (false negatives against false positives) rather than hits vs. false alarms.
using Winston ## or perhaps another plotting package
plot(r)
## The "Detection Error Tradeoff" plot, this should give a more/less straight line
detplot(r)
## compute the Area Under the ROC, should be close to 0.078
auc(r)
## define a decision cost function by its parameter p_tar=0.01, Cfa=1, Cmiss=10 (NIST SRE 2008 setting)
d = DCF(0.01, 1, 10)
## `actual costs' using a threshold of scores at -plo(d) (minus prior log odds)
plo(d)
dcf(tar, non, d=d)
## Or use a different threshold, e.g., zero
dcf(tar, non, d=d, thres=0)
## `minimal costs' using an optimal threshold
mindcf(r, d=d)
## define an array of DCFs, and compute the decision costs for these, using a threshold at -plo
d = DCF([0.001, 0.01, 0.1, 0.5, 0.9, 0.99, 0.999], 1, 1)
dcf(tar, non, d=d)
## same, but normalized to costs of decisions based on the prior alone
dcf(tar, non, d=d, norm=true)
## prior log odds, thre crucial combination of cost parameters, log(p_tar / (1-p_tar) Cmiss / Cfa)
plo(d)
## now scan the Bayes error rate (similar to the dcf above) for a range of prior log odds, and plot
## This is known as the Applied Probability of Error plot
apeplot(r)
## The area under the red curve (actual error rates), the cost of the log-likelihood ratio
cllr(tar, non)
## The area under the green curve (minimum errors), the cost of the optimal log-likelihood-ratio
minclrr(tar, non)
## Similar to APE, but normalized---a Normalized Bayes Error plot
nbeplot(r)
## Make an `LLR' plot: score-to-optimal-LLR mapping, r.θ, vs. r.llr
llrplot(r)
```

Receiver Operating Characteristic
----
A binary classifier maps an input `x` to one of two classes, `A` and `B`.  Internally, every classifier ends up producing some form of a scalar metric `s`, which can be thresholded to produce a decision.

There are two possible "senses" of this internal scalar:
- higher values of `s` indicate higher probability of `A`
- higher values of `s` indicate higher probability of `B`

There are of course also many different interpretations of the classes `A` and `B`.  For instance, in biometrics `B` could mean "same individual" and `A` "different individual".  The corresponding senses of `s` then have an interpretation
- score-like: a higher value means a better match
- distance-like: a higher value means a larger difference.

Because we want to focus in this package on a probabilistic interpretation of the scalar `s`, we take the "score-like" interpretation of `s`, i.e., higher values of `s` correspond to a higher likelihood of the class-of-interest to be associated to the input of the classifier.  If your classifier is, in fact, a distance metric `d`, you could work with `s=-d` or `s=1/d` or any other strictly decreasing function.  Alternatively, you can swap around the label of the class of interest. 

As the name suggests, a classifier is supposed to make decisions.  Decisions can be thesholded against a fixed `θ` such that:
- if `s>θ`, decide class `B`
- if `s<θ`, decide class `A`

For evaluating the classifier, we need a set of supervised trials, i.e., for each scalar score `s` we need a label indicating the true class of the trial that led to score `s`.  Because there are two classes, two types of errors can be made:
- false positives: `s>θ` while in fact the true label is `A`
- false negatives: `s<θ` while in fact the true label is `B`.

The _Receiver Operating Characteristic_ (ROC) is a graph that shows how the fractions of the false positives and false negatives change with varying `θ`, for a fixed set of scores `s`.  The structure of type `Roc` captures the essential information in a pre-processed way such that other quantities can be derived efficiently. 

Because we come from automatic speaker recognition, we tend to use the following terminology for the classes
- *target*, the higher scores, a.k.a. same source, true client, ...
- *non-target*, the lower scores, a.k.a. different source, impostor, ...

Error rates
------
In this package, we focus at analysing everything in terms of _error rates_.  Traditionally, researchers have used ROC curves where one of the axes is used to describe the complement of an error rate.  Specifically, one often sees a "true positive rate" versus a "false positive rate", where the true positive rate is the complement of the false negative rate.  There is no real objective way to judge whether one analysis is any better than the other, usually the choice is largely dependent on traditions in the area of the research you are in.

There are many different names of the error rates in different scientific disciplines.  Because I come from the area of automatic speaker recognition, the current terminology is
 - *Probability of False Alarm*: a.k.a. false alarm rate, false positive rate, false accept rate, false match rate, Type II error, 1 - true negative rate, 1- specificity
 - *Probability of Miss:*: a.k.a. miss rate, false negative rate, false reject rate, false non-match rate, Type I error, 1 - true positive rate, 1 - sensitivity, 1 - recall, 1 - verification rate, 1 - hit rate, 1 - genuine acceptance rate

We foresee that the naming of things becomes a bit more flexible in future releases of this package. 

DET plot
-------

A _detection error trade-off_ plot (Martin, 1997) is exactly the same as a ROC plot in the error domain (i.e., miss rate vs false alarm rate), but the axes have been warped according to `Φ⁻¹(x)`, the inverse of the cumulative normal distribution.  In R, this function is known as `qnorm(x)`, in Julia base this is ` √2 * erfinv(2x-1)`.  This type of plot has interesting properties
- If the distributions of target and non-target scores are both Normal, then the DET-curve is a straight line.  In practice, many detection problems give rise to more-or-less straight DET curves, and this suggests that there exists a strictly increasing warping function that can make the score distributions (more) Normal. 
- Towards better performance (lower error rates), the resolution of the graph is higher.  This makes it more easy to have multiple systems / performance characteristics over a smaller or wider range of performance in the same graph, and still be able to tell these apart.
- Conventionally, the ranges of the axes are chosen 0.1%--50%.  This makes it possible to immediately assess the overall performance based on the absolute position of the line in the graph if you have seen more DET plots in your life.
- The slope of the (straight) line corresponds to the ratio of the `σ` parameters of the underlying Normal score distributions, namely that of the non-target scores divided by that of the target scores.  Often, highly discriminative classifiers show very _flat_ curves, indicating that that target scores have a much larger variance than the non-target scores.  
- The origin of this type of plot lies in psychophysics, where graph paper with lines according to this warping was referred to as _double probability paper_.  The diagonal `y=x` in a DET plot corresponds linearly to a quantity known as `d'` (d-prime) from psychophysics, ranging from 0 at 50% error to about 6 at 0.1% error. 

### Discrete and continuous scores

There is an essential difference between discrete score (classes) and continuous scores.  For the former, trials with the same scores must be grouped before the probabilities of false alarm and miss are computed.  This results in ROC and DET plots that can have line elements that are not solely horizontal or vertical.  This is contrary to the latter case if we assume that no two scores are (coincidentally)  the same, which leads to only horizontal and vertical line segments.  This `ROCAnalysis` package makes sure that the occurrence of identical scores is treated correctly by sorting target scores before identical non-target scores, and by treating trials with scores _at_ the threshold always as errors. 

### Plot optimisation

For larget trial sets, it is very likely that in the extrems of the score distributions there is very little overlap.  This wil results in many consecutive horizontal or vertical line segments in the plot.   This `ROCAnalysis` package integrates these consecutive line segments and replaces them by a single segment, which leads to a strong reduction in complexity in further calculations and plotting.  

## Single-numbered metrics

The ROC and DET plots shows the discrimination capability of the detector as a graph.  Often one wants to summarize the plot in a single metric.  There are many ways to do this, we list some here
- *Equal Error Rate*.  This is the point in the DET or ROC where the curve crosses the `y=x` line, i.e., where the error rates are the same.  A lower EER means a better discriminating classifier.  It samples the ROC on only a single operating point, and moreover, this is an "after-the-fact" point.  This metric is insensitive to calibration, i.e., any strictly increasing function can be applied to the scores and an identical EER will be computed.  For small number of trials it makes a different how to comput the EER, often in literature this is not specified.  In this package there are the functions:
   - `eer()`: take the points where the difference between miss and false alarm rates changes sign, and interpolate to find the crossing with the `y=x` diagonal.
   - `eerch()`: compute the convex hull of the ROC, and compute the point where this crosses the `y=x`line.  This has the interpretation of _the maximum over priors of the minimum cost_, and is useful for cost function analysis.  It tends to be more stable than `eer()`.
- *Area Under the Curve*. `auc()` integrates the _Area Under the Curve_ of the ROC.  Contrary to the EER, this metric is sensitive to the entire range of operating points, but, like the EER, it is insensitive to calibration.  It can be interpreted as the probability that a random non-target score is higher than a random target score, and lower AUC indicates better discrimination.  Please note that the complement (area under the hit-rate-vs-false-alarm-rate curve) is known under the same name in other disciplines.  We believe in errors, and minimizing them.
- *Dectection Cost Function*. `dcf()`  A Detection cost function is a weighted sum of false alarms and misses.  The weights consists of separate costs for false alarms and misses, and a prior for non-targets and targets.
- *Minimum Detection Cost*. `mindcf()`  This is the minimum detection cost `dcf()` that can be obtained by varying the threshold.  It is similar to the EER in the sense that it requires setting a threshold "after the fact".  The minimum DCF is related to the points on the convex hull of the ROC curve.  
- *Cost of the Log-Likelihood-Ratio*. `cllr()` computes a normalized form of the cross-entropy between the "true posterior" (`1` for target trials and `0` for non-target trials) and the posterior resulting from the classifier score when interpreted as a likelihood ratio, and using a prior for the classes of 0.5.  This measure is _calibration sensitive_, i.e., it penalizes under- or over-confident likelihood ratios.  The minimum value is determined by the discriminative properties of the classifier, and this minimum approaches 0 for a classifier that completely separates the two classes.  A value of 1 indicates that the classifier gives no information, i.e., decisions can just as well be made based on the prior only.  A value larger than 1 indicates that making Bayes's decisions based on the classifiers score gives a higher expected cost than basing decisions on the prior alone.
- *Minimum Cllr", `mincllr()` computes the minimum attainable Cllr by warping the scores to log-likelihood-ratios while maintaining the order of the scores.  This is equivalent to determining a minimum cost for all cost functions that can be written as a linear combination of actual miss- and false-alarm-rates, and integrating these costs over cost function parameters.

Types
----
We have defined the following types:
- `TNT(tar, non)`,  A container for arrays of target and non-target scores.   Anywhere it says `tar, non` you should be able to substitute a variable of type `TNT`.
- `DCF(ptar,cfa,cmiss)`.  A container for scalars or arrays of the cost function parameters `ptar` (the prior probability of a target), `cfa` (the cost of a false alarm) and `cmiss` (the cost of a miss).  In any of the parameters are arrays, they must be of compatible size.
- `Roc(pfa, pmiss, θ, chull, llr)`.  The basic structure for storing ROC results. It contains the x- and y-coordiantes of the (error-based) ROC, together with the threshold.  Many of the "minimum" performance measures can be computed from this.  The fields in `Roc` have the following interpretation:
 - `pfa`.  The probability of a false alarm, the fraction of non-target scores _above_ (or equal to) the threshold `θ`.  This quantity is also known as false alarm/accept/positive/match rate or type II error, and `1-pfa` is known as the true negative/reject rate or specificity.
 - `pmiss`.  The probability of a miss, the fraction of target scores _below_ (or equal to) the threshold `θ`.  This quantity is also known, a.o. miss rate, false negative/reject/non-match rate, type I error, and `1-pmiss` is also known as true positive rate, sensitivity, recall, verification/hit rate, genuine acceptance rate.
 - `θ`. The threshold associated with the elements in `pfa` and `pmiss`.  It can be seen as lying between _this_ `pfa/miss` and _the next_.  Hence, there is one value less than in the `pmiss` and `pfa` arrays.
 - `chull`.  A boolean indicating whether or not the point `(pfa,pmiss))` lies on the _convex hull_.
 - `llr`.  The optimal log-likelihood-ratio associated with scores at `θ`.  


Notes
-----
This is very much work in progress.  If you stumble upon this, please drop me a line. 
