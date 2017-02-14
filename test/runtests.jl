using RSDSGE
using JLD
using Base.Test

# Run the tests of the model constructor
@testset "Model setup" begin
    include("modelsetup.jl")
end

# Check model utils work properly
@testset "Solving" begin
    include("solve.jl")
end

