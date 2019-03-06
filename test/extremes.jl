## extremes.jl  Test some extreme cases
## (c) 2014--2015 David A. van Leeuwen
##
## Licensed under the MIT software license, see LICENSE.md

## check some pathological cases where sorting and comparison of equal values plays a role

using Test

function test_both(tar, non, e1, e2)
    r = roc(tar, non)
    eer1 = eerch(r)
    eer2 = eer(tar, non)
    #println("EERch = ", eer1)
    #println("EER   = ", eer2)
    @test e1 ≈ eer1
    @test e2 ≈ eer2
end

test_both([0;2], [1], 1/3, 1/2)
test_both([1], [0;2], 1/3, 1/2)
test_both([0], [1], 1/2, 1)
test_both([1], [0], 0, 0)

test_both(collect(0:2), [1], 2/5, 1/2)
test_both(collect(0.:2), [0.5], 1/4, 1/3)
test_both(collect(0.:2), [1.5], 2/5, 2/3)

test_both([1], collect(0:2), 2/5, 1/2)
test_both([0.5], collect(0.:2), 2/5, 2/3)
test_both([1.5], collect(0.:2), 1/4, 1/3)
