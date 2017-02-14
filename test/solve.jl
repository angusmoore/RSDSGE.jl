function solve2eqNK()
    phi_H = 0.9
    phi_L = 1.5
    rho = 0.9

    parameters = ["rho","phi"]
    parameter_values = [rho, (phi_L, phi_H)]
    regime_transmatrix = [0.3 0.7; 0.3 0.7]
    vars = ["c","inflation"]
    shocks = ["e_c"]
    expressions = ["1/exp(c) = 1/exp(c(+1))*((exp(inflation))^(phi_R))/exp(inflation(+1))",
    "c = rho*c(-1) + e_c"]

    twoEQNK = RSDSGEModel(vars,shocks,parameters,expressions,parameter_values,regime_transmatrix)
    return solve(twoEQNK)
end

function solveRBCfromFRWZ()
    beta = 0.9976
    mu = (0.005,0.00004)
    sigma = 0.0002
    alpha = 0.33
    delta = 0.025
    transmatrix = [0.9 0.1; 0.1 0.9]

    vars = ["c","k","eps"]
    parameters = ["betta","mu","sigma","alpha","delta"]
    parameter_values = [beta,mu,sigma,alpha,delta]

    shocks = ["e_eps"]

    expressions = ["1/c = betta*1/c(+1)*exp((mu_R + sigma*eps)/(alpha-1))*(alpha*exp(mu_RP + sigma*eps(+1))*k^(alpha-1)+1-delta)",
    "c + k*exp((mu_R + sigma*eps)/(1-alpha)) = exp(mu_R+sigma*eps)*(k(-1))^(alpha) + (1-delta)*k(-1)",
    "eps = e_eps"]

    FRWZ_RBC = RSDSGEModel(vars,shocks,parameters,expressions,parameter_values,transmatrix,ones(Float64,3))
    return solve(FRWZ_RBC)
end

NKsol = solve2eqNK()
FRWZsol = solveRBCfromFRWZ()

@load "modelsolutions.jld"

@testset "Number of solutions" begin
    @test length(NKsol)==1
    @test length(FRWZsol) == 2
end

@testset "MX matrices" begin
    @test_approx_eq NKsol[1].MX NKsolMX
    @test_approx_eq FRWZsol[1].MX FRWZsolMX1
    @test_approx_eq FRWZsol[2].MX FRWZsolMX2
end

@testset "ME matrices" begin
    @test_approx_eq NKsol[1].ME NKsolME
    @test_approx_eq FRWZsol[1].ME FRWZsolME1
    @test_approx_eq FRWZsol[2].ME FRWZsolME2
end

@testset "MC matrices" begin
    @test_approx_eq NKsol[1].MC zeros(Float64,2,1,2)
    @test_approx_eq FRWZsol[1].MC FRWZsolMC1
    @test_approx_eq FRWZsol[2].MC FRWZsolMC2
end


