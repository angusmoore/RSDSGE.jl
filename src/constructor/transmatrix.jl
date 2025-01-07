function getergodic(transmatrix::Array{Float64,2})
    D = LinearAlgebra.eigvals(transmatrix')
    V = LinearAlgebra.eigvecs(transmatrix')
    _, eigindex = findmin(abs.(D .- 1.0))
    ergodic = vec(V[:,eigindex]/sum(V[:,eigindex]))
    if any(abs.(imag(ergodic)).>1e-16)
        error("Imaginary ergodic distribution")
    else
        return real(ergodic)
    end
end

function createtransmatrix_sym(numregimes::Int)
    return Symbolics.Num[Symbolics.variable("TRANSITIONPROBABILITY_RX_RP$(j)") for j in 1:numregimes]
end

function TransitionMatrix(trans::Array{Float64,2})
    # Sanity check the regime transition matrix
    numregimes = size(trans,1)
    if size(trans,2)!=size(trans,1)
        error("Regime transition matrix is not square.")
    end

    if any((abs.(sum(trans, dims=2)) .- 1.0) .> 1e-12)
        error("Rows of the regime transition matrix do not sum to 1.")
    end
    # Create a sympy version of the regime trans matrix
    syms = createtransmatrix_sym(numregimes)

    # Get the ergodic distribution of the markovian regime process
    ergodic = getergodic(trans)
    return TransitionMatrix(trans,syms,ergodic,Symbolics.variable("TRANSITIONPROBABILITY"))
end
