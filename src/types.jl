type RegimeParameterDecomposition
    bar::Float64
    hat::Array{Float64,1}
end

type Block
    size::Int
    equations::Array{Bool,1}
    vars::Array{Bool,1}
    states::Array{Bool,1}
end

immutable Meta
    numstates::Int
    numvars::Int
    numshocks::Int
    numregimes::Int
    numparameters::Int
    numstatic::Int
    numexo::Int
end

immutable Parameters
    names::Array{String,1}
    values::Array{Float64,2}
    dictionary::Dict{String,Int}
    isswitching::Array{Bool,1}
    affectsSS::Array{Bool,1}
    decomposition::Array{RegimeParameterDecomposition,1}
    values_withbar::Array{Float64,2}
    generic_sym::Array{SymPy.Sym,2}
    generic_hats::Array{SymPy.Sym,2}
    bar_sym::Array{SymPy.Sym,1}
    specificregime_sym::Array{SymPy.Sym,2}
end

immutable Variables
    names::Array{String,1}
    dictionary::Dict{String,Int}
    isstate::Array{Bool,1}
    isstatic::Array{Bool,1}
    isexo::Array{Bool,1}
    states_indices::Array{Int,1}
    states_sym::Array{SymPy.Sym,1}
    leads_sym::Array{SymPy.Sym,1}
    contemps_sym::Array{SymPy.Sym,1}
end

immutable Shocks
    names::Array{String,1}
    dictionary::Dict{String,Int}
    syms::Array{SymPy.Sym,1}
end

immutable Equations
    stringversions::Array{String,1}
    isstatic::Array{Bool,1}
    isexo::Array{Bool,1}
    syms::Array{SymPy.Sym,1}
end

type SteadyState # not immutable, because I need to wholesale reassign the steady state values all the time
    system::Array{Function,1}
    values::Array{Float64,1}
end

immutable TransitionMatrix
    values::Array{Float64,2}
    genericsyms::Array{SymPy.Sym,1}
    ergodic::Array{Float64,1}
    completelygeneric::SymPy.Sym
end

immutable PerturbationSystem
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
    isJ::Array{Bool,1}
    isY::Array{Bool,1}
    isW::Array{Bool,1}
    isZ::Array{Bool,1}
    isX::Array{Bool,1}
    blocks::Array{Block,1}
end

type RSDSGEModel
    meta::Meta
    parameters::Parameters
    vars::Variables
    shocks::Shocks
    equations::Equations
    steadystate::SteadyState
    transmatrix::TransitionMatrix
    perturbationsystem::PerturbationSystem
end