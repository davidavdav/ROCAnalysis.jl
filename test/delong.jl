@testset "test delong" begin
    target = [0,0,1,1,1]
    modela = [0.1, 0.2, 0.6, 0.7, 0.8]
    modelb = [0.3,0.6,0.2,0.7,0.9]
    nontara, tara = modela[1:2], modela[3:end]
    nontarb, tarb = modelb[1:2], modelb[3:end]

    @test round(ROCAnalysis.delong_test(tara, nontara, tarb, nontarb), digits=4) == 0.3173 

    target = [0,0,0,0,0,0,1,1,1,1,1,1,1]
    modela = [0.1,0.2,0.05,0.3,0.1,0.6,0.6,0.7,0.8,0.99,0.8,0.67,0.5]
    modelb = [0.3,0.6,0.2,0.1,0.1,0.9,0.23,0.7,0.9,0.4,0.77,0.3,0.89]
    nontara, tara = modela[1:6], modela[7:end]
    nontarb, tarb = modelb[1:6], modelb[7:end]

    @test round(ROCAnalysis.delong_test(tara, nontara, tarb, nontarb), digits=5) == 0.09453
end
