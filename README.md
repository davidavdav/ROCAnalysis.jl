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

Synopsis
----
```julia
## Produce some well-calibrated log-likelihood-ratio scores:
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

Receiver Operating Cheracteristic
----
A binary classifier maps an input `x` to one of two classes, `A` and `B`.  Internally, every classifier ends up producing some form of a scalar metric `s`, which can be thresholded to produce a decision.

There are two possible "senses" of this internal scalar:
- higher values of `s` indicate higher probability of `A`
- higher values of `s` indicate higher probability of `B`
There are of course also many different interpretations of the classes `A` and `B`.  For instance, in biometrics `A` could mean "same indicidual" and `B` "different individual".  The corresponding senses of `s` then have an interpretation
- score-like: a higher value means a better match
- distance-like: a higher value means a larger difference.

Because of this package we want to focus on a probabilistic interpretation of the scalar `s`, we take the "score-like" interpretation of `s`, i.e., higher values of `s` correspond to a higher likelihood of the class-of-interest to be associated to the input of the classifier.  If your classifier is, in fact, a distance metric `d`, you could work with `s=-d` or `s=1/d` or any other strictly decreasing function.



Error rates
------

Notes
-----
This is very much work in progress.  If you stumble upon this, please drop me a line. 
