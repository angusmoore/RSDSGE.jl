using RSDSGE
using JLD2
using Test

using Aqua
Aqua.test_all(RSDSGE)

# Run the tests of the model constructor
@testset "Model setup" begin
    include("test-modelsetup.jl")
end

# Check model solving works properly
@testset "Solving" begin
    include("test-solve.jl")
end

# Check model utils
@testset "Model utils" begin
    include("test-modelutils.jl")
end
