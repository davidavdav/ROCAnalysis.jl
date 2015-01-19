## dcf.jl  Tests related to decision cost functions

using Base.Test

x = readtable("ru.2009.table.gz", separator='\t')
tnt = TNT(x)
tar, non = tnt.tar, tnt.non
ruc = roc(tnt, collapse=false)
r = roc(tnt)

## scalar
d = DCF(0.01, 1, 10)            # nist 2008
@test_approx_eq oeff(d) 1/9.9
@test_approx_eq peff(d) 1/10.9
@test_approx_eq plo(d) -log(9.9)

lo = plo(d)
b = 0.03291547867387523
@test_approx_eq ber(x,lo) b
@test_approx_eq ber(tnt,lo) b
@test_approx_eq ber(tar,non,lo) b
@test_approx_eq ber(ruc,lo) b

m = 0.03173675448189074
@test_approx_eq minber(x,lo) m
@test_approx_eq minber(tnt,lo) m
@test_approx_eq minber(tar,non,lo) m
@test_approx_eq minber(r,lo) m

for norm in (false, true)
    c = 0.035877871754524004 * 10^norm
    @test_approx_eq dcf(x,d,norm=norm) c
    @test_approx_eq dcf(tnt,d,norm=norm) c
    @test_approx_eq dcf(tar,non,d,norm=norm) c
    @test_approx_eq dcf(ruc,d,norm=norm) c

    m = 0.034593062385260914 * 10^norm
    @test_approx_eq mindcf(x,d,norm=norm) m
    @test_approx_eq mindcf(tnt,d,norm=norm) m
    @test_approx_eq mindcf(tar,non,d,norm=norm) m
    @test_approx_eq mindcf(r,d,norm=norm) m
end

## array
lo = [-7:0.01:7]
dc = DCF(sigmoid(lo), 1, 1)

@test_approx_eq_eps ber(x,lo) ber(ruc,lo) 1e-5
@test ber(tnt,lo) == ber(tar,non,lo)

@test minber(x,lo) == minber(tnt,lo) == minber(tar,non,lo) == minber(r,lo)
@test all(minber(r,lo) .≤ ber(ruc,lo))

for norm in (false, true)
    @test_approx_eq dcf(x,dc,norm=norm) dcf(ruc,dc,norm=norm)
    @test dcf(x, dc, norm=norm) == dcf(tnt, dc, norm=norm) == dcf(tar, non, dc, norm=norm)
    @test mindcf(x, dc, norm=norm) == mindcf(tnt, dc, norm=norm) == mindcf(tar, non, dc, norm=norm) == mindcf(r, dc, norm=norm)
    @test all(mindcf(r, dc, norm=norm) .≤ dcf(ruc, dc, norm=norm))
end
