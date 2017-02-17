path=dirname(@__FILE__)
f = load("$path/modelsolutions.jld")

NKsol=solve(twoEQNK)
FRWZsol=solve(FRWZ)

@testset "Number of solutions" begin
    @test length(NKsol)==1
    @test length(FRWZsol) == 2
end

@testset "MX matrices" begin
    @test isapprox(NKsol[1].MX, f["NKsolMX"])
    @test isapprox(FRWZsol[1].MX, f["FRWZsolMX1"])
    @test isapprox(FRWZsol[2].MX, f["FRWZsolMX2"])
end

@testset "ME matrices" begin
    @test isapprox(NKsol[1].ME, f["NKsolME"])
    @test isapprox(FRWZsol[1].ME, f["FRWZsolME1"])
    @test isapprox(FRWZsol[2].ME, f["FRWZsolME2"])
end

@testset "MC matrices" begin
    @test isapprox(NKsol[1].MC, zeros(Float64,2,1,2))
    @test isapprox(FRWZsol[1].MC, f["FRWZsolMC1"])
    @test isapprox(FRWZsol[2].MC, f["FRWZsolMC2"])
end


