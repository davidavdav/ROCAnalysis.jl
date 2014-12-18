ROC.jl
======
[![Build Status](https://travis-ci.org/davidavdav/ROC.jl.svg?branch=master)](https://travis-ci.org/davidavdav/ROC.jl)


Receiver Operating Characteristics and functions for evaluation probabilistic binary classifiers. 

Please note there is an alternative implementation under [the same name](https://github.com/diegozea/ROC.jl).  This implementation is more geared towards: 
 - large amounts of data, with efficient ROC statistics calculation
 - Detection Error Trade-off (DET) analysis
 - ROC convex hull computation and EER-interpretation 
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

Notes
-----
This is very much work in progress.  If you stumble upon this, please drop me a line. 
