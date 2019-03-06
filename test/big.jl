## big.jl  A test with a larger set of data
## (c) 2014--2015 David A. van Leeuwen
##
## Licensed under the MIT software license, see LICENSE.md

using Test
import Random

## Note that due to the random sample, some of these test may fail every once in a while,
## So let's set a seed
Random.seed!(1)

## test some statistics on large data sets
## calibrated scores (i.e., these scores can be interpreted as log-likelihood-ratios)
tar =  2 .+ 2randn(1_000_000)
non = -2 .+ 2randn(1_000_000)

## eer
@test eer(tar, non) ≈ pnorm(-1) atol=0.01

## eerch
@time r = roc(tar, non)
@test eerch(r) ≈ pnorm(-1) atol=0.01

## AUC
pauc = 0.0786496
@test auc(r) ≈ pauc atol=0.01
auc(r, pfa=0.1)
AUC(r, pfa=0.1)

## Cllr
cllr2 = 0.5140558 ## Cllr for distributions with μ = ±2 and σ = 2, i.e., d'=2
@test cllr(tar, non) ≈ cllr2 atol=0.01
@test mincllr(tar, non) ≈ cllr2 atol=0.01
