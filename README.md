ROC.jl
======
[![Build Status](https://travis-ci.org/davidavdav/ROC.jl.svg?branch=master)](https://travis-ci.org/davidavdav/ROC.jl)


Receiver Operating Characteristics and functions for evaluation probabilistic binary classifiers. 

Please note there is an alternative implementation under [the same name](https://github.com/diegozea/ROC.jl).

Our implementation is more geared towards:
 - large amounts of data, with efficient ROC statistics calculation
 - Cost function analysis
 - Detection Error Trade-off (DET) analysis
 - ROC convex hull computation, analysis and EER-interpretation
 - Optimal Likelihood Ratio computation

The development roadmap is largely based on the functionality in a similar R package [ROC](https://github.com/davidavdav/ROC). 

Synopsis
----
```julia
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
plot(r)
## this should give a more/less straight line
detplot(r)
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

The _Receiver Operating Characteristic_ (ROC) is a graph that shows how the fractions of the false positives and false negatives change with varying `θ`, for a fixed set of scores `s`.  

Error rates
------
In this package, we focus at analysing everything in terms of _error rates_.  Traditionally, researchers have used ROC curves where one of the axes is used to describe the complement of an error rate.  Specifically, one often sees a "true positive rate" versus a "false positive rate", where the true positive rate is the complement of the false negative rate.  There is no real objective way to judge whether one analysis is any better than the other, usually the choice is largely dependent on traditions in the area of the research you're in.

There are many different names of the error rates in different scientific disciplines.  Because I come from the area of automatic speaker recognition, the current terminology is
 - *Probability of False Alarm*: false alarm rate, false positive rate, 1 - true negative rate
 - *Probability of Miss: miss rate*, false negative rate, 1 - hit rate, 1 - true positive rate

Notes
-----
This is very much work in progress.  If you stumble upon this, please drop me a line. 
