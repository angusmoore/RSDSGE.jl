module RSDSGE

using Symbolics: Symbolics, parse_expr_to_symbolic
using Groebner: Groebner
import Base: copy
using NLsolve: NLsolve
using SparseArrays: SparseArrays, SparseMatrixCSC
using LinearAlgebra: LinearAlgebra
using Polynomials: Polynomials
using SymbolicUtils: SymbolicUtils

baremodule EqEvalModule
using Symbolics: Symbolics
using Base: *, +, -, /, ^, !, ==, exp, log
import Base.@eval
import Base.Iterators
import Base.any
import Base.enumerate
import Base.map

end

include("constructor/modeltypes.jl")
include("constructor/utils.jl")
include("constructor/constructor.jl")
include("constructor/parsefile.jl")
include("constructor/perturbation.jl")
include("modelutils.jl")
include("solve/solve.jl")

# Functions to be exported
export RSDSGEModel, RSDSGESolution, solve, updateparameters!, findSS!

end
