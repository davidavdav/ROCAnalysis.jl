## compare my chull implementation with that of the general qhull library

function testqhull(tar, non)
    xo, tc = sortscores(tar, non)
    θ, tc, nc, = binscores(xo, tc)
    pfa = [1, 1 - cumsum(nc)//length(non)]
    pmiss = [0, cumsum(tc)//length(tar)]
    @time ch, llr = chllr(tc, nc, xo, laplace=false, fast=false)
    @time myqh = rochull(pfa, pmiss)
    pfa, pmiss, find(ch), myqh
end

function test()
    ntar = rand(1:10000)
    nnon = rand(1:10000)
    tar = 1 + randn(ntar)
    non = -1 + randn(nnon)
    pfa, pmiss, a, b = testqhull(tar, non)
    if a != b
        println("Difference found:")
        if length(a) > length(b)
            println("qhull more than rochull")
            println(setdiff(a,b))
        else
            println("rochull more than qhull")
            println(setdiff(b,a))
        end
    end
    return pfa, pmiss, a, b
end

function prof()
    N = 1000_000
    tar = 2 + 2randn(N)
    non = -2 + 2randn(N)
    Profile.clear()
    @profile roc(tar, non)
    Profile.print(format=:flat)
end
