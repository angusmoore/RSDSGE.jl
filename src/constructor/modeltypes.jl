mutable struct RegimeParameterDecomposition
    bar::Float64
    hat::Array{Float64,1}
end

struct Block
    size::Int
    equations::BitArray{1}
    vars::BitArray{1}
    states::BitArray{1}
end

struct Meta
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
    isswitching::BitArray{1}
    affectsSS::BitArray{1}
    decomposition::Array{RegimeParameterDecomposition,1}
    values_withbar::Array{Float64,2}
    generic_sym::Array{SymPy.Sym,2}
    generic_hats::Array{SymPy.Sym,2}
    bar_sym::Array{SymPy.Sym,1}
    specificregime_sym::Array{SymPy.Sym,2}
end

struct Variables
    names::Array{String,1}
    dictionary::Dict{String,Int}
    isstate::BitArray
    isstatic::BitArray
    isexo::BitArray
    states_indices::Array{Int,1}
    states_sym::Array{SymPy.Sym,1}
    leads_sym::Array{SymPy.Sym,1}
    contemps_sym::Array{SymPy.Sym,1}
end

struct Shocks
    names::Array{String,1}
    dictionary::Dict{String,Int}
    syms::Array{SymPy.Sym,1}
end

struct Equations
    stringversions::Array{String,1}
    isstatic::BitArray{1}
    isexo::BitArray{1}
    syms::Array{SymPy.Sym,1}
end

mutable struct SteadyState # mutable, because I need to wholesale reassign the steady state values all the time
    system::Array{Function,1}
    values::Array{Float64,1}
end

struct TransitionMatrix
    values::Array{Float64,2}
    genericsyms::Array{SymPy.Sym,1}
    ergodic::Array{Float64,1}
    completelygeneric::SymPy.Sym
end

struct PerturbationSystem
    A1::SparseMatrixCSC{Function}
    A2::SparseMatrixCSC{Function}
    A3::SparseMatrixCSC{Function}
    A4::SparseMatrixCSC{Function}
    B1::SparseMatrixCSC{Function}
    B2::SparseMatrixCSC{Function}
    B3::SparseMatrixCSC{Function}
    B4::SparseMatrixCSC{Function}
    B5::SparseMatrixCSC{Function}
    C1_states::SparseMatrixCSC{Function}
    C2_states::SparseMatrixCSC{Function}
    C_shocks::SparseMatrixCSC{Function}
    C_constant::SparseMatrixCSC{Function}
    C_constant_p::SparseMatrixCSC{Function}
    Xs::Array{SymPy.Sym,3}
    Xx::Array{SymPy.Sym,2}
    Xz::Array{SymPy.Sym,2}
    Js::Array{SymPy.Sym,3}
    Jx::Array{SymPy.Sym,2}
    Jz::Array{SymPy.Sym,2}
    isJ::BitArray{1}
    isY::BitArray{1}
    isW::BitArray{1}
    isZ::BitArray{1}
    isX::BitArray{1}
    blocks::BitArray{1}
end

struct RSDSGEModel
    meta::Meta
    parameters::Parameters
    vars::Variables
    shocks::Shocks
    equations::Equations
    steadystate::SteadyState
    transmatrix::TransitionMatrix
    perturbationsystem::PerturbationSystem
end
