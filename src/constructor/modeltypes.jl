mutable struct RegimeParameterDecomposition
    bar::Float64
    hat::Array{Float64,1}
end

struct Block
    size::Int
    equations::Vector{Bool}
    vars::Vector{Bool}
    states::Vector{Bool}
end

struct ModelMetaInfo
    numstates::Int
    numvars::Int
    numshocks::Int
    numregimes::Int
    numparameters::Int
    numstatic::Int
    numexo::Int
end

struct Parameters
    names::Array{String,1}
    values::Array{Float64,2}
    dictionary::Dict{String,Int}
    isswitching::Vector{Bool}
    affectsSS::Vector{Bool}
    decomposition::Array{RegimeParameterDecomposition,1}
    values_withbar::Array{Float64,2}
    generic_sym::Array{Symbolics.Num,2}
    generic_hats::Array{Symbolics.Num,2}
    bar_sym::Array{Symbolics.Num,1}
    specificregime_sym::Array{Symbolics.Num,2}
end

struct Variables
    names::Array{String,1}
    dictionary::Dict{String,Int}
    isstate::BitArray
    isstatic::BitArray
    isexo::BitArray
    states_indices::Array{Int,1}
    states_sym::Array{Symbolics.Num,1}
    leads_sym::Array{Symbolics.Num,1}
    contemps_sym::Array{Symbolics.Num,1}
end

struct Shocks
    names::Array{String,1}
    dictionary::Dict{String,Int}
    syms::Array{Symbolics.Num,1}
end

struct Equations
    stringversions::Array{String,1}
    isstatic::Vector{Bool}
    isexo::Vector{Bool}
    syms::Array{Symbolics.Num,1}
end

mutable struct SteadyState # mutable, because I need to wholesale reassign the steady state values all the time
    system::Array{Function,1}
    values::Array{Float64,1}
end

struct TransitionMatrix
    values::Array{Float64,2}
    genericsyms::Array{Symbolics.Num,1}
    ergodic::Array{Float64,1}
    completelygeneric::Symbolics.Num
end

struct PerturbationSystem
    A1::Matrix{Function}
    A2::Matrix{Function}
    A3::Matrix{Function}
    A4::Matrix{Function}
    B1::Matrix{Function}
    B2::Matrix{Function}
    B3::Matrix{Function}
    B4::Matrix{Function}
    B5::Matrix{Function}
    C1_states::Matrix{Function}
    C2_states::Matrix{Function}
    C_shocks::Matrix{Function}
    C_constant::Matrix{Function}
    C_constant_p::Matrix{Function}
    Xs::Array{Symbolics.Num,3}
    Xx::Array{Symbolics.Num,2}
    Xz::Array{Symbolics.Num,2}
    Js::Array{Symbolics.Num,3}
    Jx::Array{Symbolics.Num,2}
    Jz::Array{Symbolics.Num,2}
    isJ::Vector{Bool}
    isY::Vector{Bool}
    isW::Vector{Bool}
    isZ::Vector{Bool}
    isX::Vector{Bool}
    blocks::Vector{Bool}
end

struct RSDSGEModel
    meta::ModelMetaInfo
    parameters::Parameters
    vars::Variables
    shocks::Shocks
    equations::Equations
    steadystate::SteadyState
    transmatrix::TransitionMatrix
    perturbationsystem::PerturbationSystem
end
