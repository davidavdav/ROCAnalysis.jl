## dcf.jl  Tests related to decision cost functions

using Test
import CSV
import GZip

x = GZip.open("ru.2009.table.gz") do fd
    CSV.read(fd, delim='\t', truestrings=["TRUE"], falsestrings=["FALSE"])
end
tnt = TNT(x)
tar, non = tnt.tar, tnt.non
ruc = roc(tnt, collapse=false)
r = roc(tnt)

## scalar
d = DCF(0.01, 1, 10)            # nist 2008
@test oeff(d) ≈ 1/9.9
@test peff(d) ≈ 1/10.9
@test plo(d) ≈ -log(9.9)
## print
show(stdout, MIME"text/plain"(), d)
println([d,d])

lo = plo(d)
b = 0.03291547867387523
@test ber(x, lo) ≈ b
@test ber(tnt, lo) ≈ b
@test ber(tar, non, lo) ≈ b
@test ber(ruc, lo) ≈ b

m1 = 0.03173675448189074
@test minber(x, lo) ≈ m1
@test minber(tnt, lo) ≈ m1
@test minber(tar, non, lo) ≈ m1
@test minber(r, lo) ≈ m1

## set a global DCF as well
setdcf(cmiss=10)                # evalita 2009

for norm in (false, true)
    c = 0.035877871754524004 * 10^norm
    @test dcf(x, d=d, norm=norm) ≈ c
    @test dcf(tnt, d=d, norm=norm) ≈ c
    @test dcf(tar, non, d=d, norm=norm) ≈ c
    @test dcf(ruc, d=d, norm=norm) ≈ c

    m2 = 0.034593062385260914 * 10^norm
    @test mindcf(x, d=d, norm=norm) ≈ m2
    @test mindcf(tnt, d=d, norm=norm) ≈ m2
    @test mindcf(tar, non, d=d, norm=norm) ≈ m2
    @test mindcf(r, d=d, norm=norm) ≈ m2

    ## global DCF
    c = 0.23094942466561763084 * 2^norm
    @test dcf(x, norm=norm) ≈ c
    @test dcf(tnt, norm=norm) ≈ c
    @test dcf(tar, non, norm=norm) ≈ c
    @test dcf(ruc, norm=norm) ≈ c

    m3 = 0.21843024685287174003 * 2^norm
    @test mindcf(x, norm=norm) ≈ m3
    @test mindcf(tnt, norm=norm) ≈ m3
    @test mindcf(tar, non,norm=norm) ≈ m3
    @test mindcf(ruc, norm=norm) ≈ m3
end

## array
lo = collect(-7:0.01:7)
dc = DCF(sigmoid.(lo), 1, 1)

@test ber(x, lo) ≈ ber(ruc, lo) atol=1e-5
@test ber(tnt, lo) == ber(tar, non, lo)

@test minber(x, lo) == minber(tnt, lo) == minber(tar, non, lo) == minber(r, lo)
@test all(minber(r, lo) .≤ ber(ruc, lo))

for norm in (false, true)
    @test dcf(x, d=dc, norm=norm) ≈ dcf(ruc, d=dc, norm=norm)
    @test dcf(x, d=dc, norm=norm) == dcf(tnt, d=dc, norm=norm) == dcf(tar, non, d=dc, norm=norm)
    @test mindcf(x, d=dc, norm=norm) == mindcf(tnt, d=dc, norm=norm) == mindcf(tar, non, d=dc, norm=norm) == mindcf(r, d=dc, norm=norm)
    @test all(mindcf(r, d=dc, norm=norm) .≤ dcf(ruc, d=dc, norm=norm))
end
