struct EvaluatedMatrices
    A1::Array{Float64,2}
    A2::Array{Float64,2}
    A3::Array{Float64,2}
    A4::Array{Float64,2}
    B1::Array{Float64,2}
    B2::Array{Float64,2}
    B3::Array{Float64,2}
    B4::Array{Float64,2}
    B5::Array{Float64,2}
    C1_states::Array{Float64,2}
    C2_states::Array{Float64,2}
    C_shocks::Array{Float64,2}
    C_constant::Array{Float64,2}
end

struct StaticSub
    YA1::Array{Float64,2}
    YA2::Array{Float64,2}
    YA3::Array{Float64,2}
    YA4::Array{Float64,2}
    YB1::Array{Float64,2}
    YB2::Array{Float64,2}
    YB3::Array{Float64,2}
    YB4::Array{Float64,2}
    YC1::Array{Float64,2}
    YC2::Array{Float64,2}
end

struct ExoSolution
    W::Array{Float64,2}
    Z::Array{Float64,2}
end

struct DynamicSolution
    X::Array{Float64,2}
    J::Array{Float64,2}
end

struct StaticSolution
    Y::Array{Float64,2}
end

struct ShocksSolution
    orderedME::Array{Float64,2}
end

struct ConstantSolution
    orderedMC::Array{Float64,2}
end

struct RSDSGESolution
    MX::Array{Float64,3}
    ME::Array{Float64,3}
    MC::Array{Float64,3}
end

function EvaluatedMatrices(model::RSDSGEModel,r::Int)
    A1 = createstacked(model,model.perturbationsystem.A1,r)
    A2 = createstacked(model,model.perturbationsystem.A2,r)
    A3 = createstacked(model,model.perturbationsystem.A3,r)
    A4 = createstacked(model,model.perturbationsystem.A4,r)
    B1 = createsummed(model,model.perturbationsystem.B1,r)
    B2 = createsummed(model,model.perturbationsystem.B2,r)
    B3 = createsummed(model,model.perturbationsystem.B3,r)
    B4 = createsummed(model,model.perturbationsystem.B4,r)
    B5 = createsummed(model,model.perturbationsystem.B5,r)
    C1 = createsummed(model,model.perturbationsystem.C1_states,r)
    C2 = createsummed(model,model.perturbationsystem.C2_states,r)
    C_shocks = createsummed(model,model.perturbationsystem.C_shocks,r)
    C_constant = createstacked(model,model.perturbationsystem.C_constant,r)
    C_constant_p = createstacked(model,model.perturbationsystem.C_constant_p,r)
    multiplybyhats!(C_constant,model,r,true)
    multiplybyhats!(C_constant_p,model,r,false)
    # Collapse the constant matrices
    C_constant += C_constant_p
    C_constant = sum(C_constant,dims=2)
    return EvaluatedMatrices(A1,A2,A3,A4,B1,B2,B3,B4,B5,C1,C2,C_shocks,C_constant)
end

function evaluatematrices(model::RSDSGEModel)
    out = Vector{EvaluatedMatrices}()
    for r in 1:model.meta.numregimes
        push!(out,EvaluatedMatrices(model,r))
    end
    return out
end

function createarglistvalues(model::RSDSGEModel,r::Integer,rp::Integer)
    arg = Vector{Float64}(undef, 1)
    arg = [model.transmatrix.values[r,rp]] # Transition prob first
    append!(arg,model.steadystate.values) # Add in the steady state values
    # Non-switching vars first
    append!(arg,model.parameters.values_withbar[r, .!model.parameters.isswitching])
    # For parameters that affect SS, just need the bar version
    append!(arg,model.parameters.values_withbar[r,model.parameters.affectsSS])
    # Finally, for parameters that switch but _don't_ affect SS, need both R and and RP version
    append!(arg,model.parameters.values_withbar[r,model.parameters.isswitching .& .!model.parameters.affectsSS])
    append!(arg,model.parameters.values_withbar[rp,model.parameters.isswitching .& .!model.parameters.affectsSS])
    return arg
end

function evaluate(model::RSDSGEModel,fmatrix::Matrix{Function},r::Integer,rp::Integer)
    args = createarglistvalues(model,r,rp)
    out = zeros(Float64,size(fmatrix)...)
    for i in eachindex(fmatrix)
        out[i] = fmatrix[i](args...)
    end
    return out
end

function createsummed(model::RSDSGEModel,fmatrix,r::Integer)
    out = zeros(Float64,size(fmatrix)...)
    if !isempty(fmatrix)
        for rp in 1:model.meta.numregimes
            out += evaluate(model,fmatrix,r,rp)
        end
    end
    return out
end

function createstacked(model::RSDSGEModel,fmatrix::Matrix{Function},r::Integer)
    if !isempty(fmatrix)
        out = zeros(size(fmatrix,1),0)
        for rp in 1:model.meta.numregimes
            out = hcat(out,evaluate(model,fmatrix,r,rp))
        end
    else
        out = zeros(size(fmatrix,1),size(fmatrix,2)*model.meta.numregimes)
    end
    return out
end

function multiplybyhats!(C_constant,model,r,today)
    # This function goes through the C_constant array and multiplies each element by it's relevant hat
    # since that is the derivative of the parameter wrt perturbationconstant (and the element already in the matrix
    # is the derivative of the system wrt the parameter
    for rp in 1:model.meta.numregimes
        i = 0
        for p in 1:model.meta.numparameters
            if model.parameters.affectsSS[p]
                i+=1
                x = (rp-1)*count(!iszero, model.parameters.affectsSS) + i
                for eq in 1:size(C_constant,1)
                    if today
                        C_constant[eq,x] = C_constant[eq,x]*model.parameters.decomposition[p].hat[r]
                    else
                        C_constant[eq,x] = C_constant[eq,x]*model.parameters.decomposition[p].hat[rp]
                    end
                end
            end
        end
    end
    return
end

function substatics(model::RSDSGEModel,m,r)
    A1 = m.A1[model.equations.isstatic,:]
    A2 = m.A2[model.equations.isstatic,:]
    A3 = m.A3[model.equations.isstatic,:]
    A4 = m.A4[model.equations.isstatic,:]
    B1 = m.B1[model.equations.isstatic,:]
    B2 = m.B2[model.equations.isstatic,:]
    B3 = m.B3[model.equations.isstatic,:]
    B4 = m.B4[model.equations.isstatic,:]
    B5 = m.B5[model.equations.isstatic,:]
    C1 = m.C1_states[model.equations.isstatic,:]
    C2 = m.C2_states[model.equations.isstatic,:]
    return StaticSub(-B5\A1,-B5\A2,-B5\A3,-B5\A4,-B5\B1,-B5\B2,-B5\B3,-B5\B4,-B5\C1,-B5\C2)
end

function substatics(model::RSDSGEModel,m)
    return StaticSub[substatics(model,m[r],r) for r in 1:model.meta.numregimes]
end

function solveexos(model::RSDSGEModel,m::EvaluatedMatrices,staticsol::StaticSub,r::Int)
    # These shouldn't be needed
    A1 = m.A1[model.equations.isexo,:]
    A2 = m.A2[model.equations.isexo,:]
    A3 = m.A3[model.equations.isexo,:]
    A4 = m.A4[model.equations.isexo,:]
    B1 = m.B1[model.equations.isexo,:]
    B2 = m.B2[model.equations.isexo,:]

    B3 = m.B3[model.equations.isexo,:]
    B4 = m.B4[model.equations.isexo,:]
    B5 = m.B5[model.equations.isexo,:]
    C1 = m.C1_states[model.equations.isexo,:]
    C2 = m.C2_states[model.equations.isexo,:]

    A1 += B5*staticsol.YA1
    A2 += B5*staticsol.YA2
    A3 += B5*staticsol.YA3
    A4 += B5*staticsol.YA4
    B1 += B5*staticsol.YB1
    B2 += B5*staticsol.YB2
    B3 += B5*staticsol.YB3
    B4 += B5*staticsol.YB4
    C1 += B5*staticsol.YC1
    C2 += B5*staticsol.YC2

    # These should all be non zero. This is just double checking that.
    @assert count(!iszero, A1) == 0 "Non-zero A1 entries in exo block"
    @assert count(!iszero, A2) == 0 "Non-zero A2 entries in exo block"
    @assert count(!iszero, A3) == 0 "Non-zero A3 entries in exo block"
    @assert count(!iszero, A4) == 0 "Non-zero A4 entries in exo block"
    @assert count(!iszero, B1) == 0 "Non-zero B1 entries in exo block"
    @assert count(!iszero, B3) == 0 "Non-zero B3 entries in exo block"

    sol = -hcat(B2,B4)\hcat(C1,C2)
    W = sol[1:count(!iszero, model.perturbationsystem.isW),:]
    Z = sol[(1+count(!iszero, model.perturbationsystem.isW)):end,:]
    return ExoSolution(W,Z)
end

function solveexos(model::RSDSGEModel,m,staticsubs)
    return ExoSolution[solveexos(model,m[r],staticsubs[r],r) for r in 1:model.meta.numregimes]
end

function arrangeexomatrices(model::RSDSGEModel,exosols)
    nz = count(!iszero, model.perturbationsystem.isZ)
    nw = count(!iszero, model.perturbationsystem.isW)
    nx = count(!iszero, model.perturbationsystem.isX)
    Wz = zeros(Float64,0,nz)
    Zz = zeros(Float64,0,nz)
    for r in 1:model.meta.numregimes
        Wz = vcat(Wz,exosols[r].W[:,(1+nx):end])
        Zz = vcat(Zz,exosols[r].Z[:,(1+nx):end])
    end
    return Wz,Zz
end

function createsystem(model::RSDSGEModel,m::EvaluatedMatrices,staticsol::StaticSub,exosol::ExoSolution,Wz::Matrix{Float64},Zz::Matrix{Float64},dynamics::Vector{Bool},r::Int)
    A1 = m.A1[dynamics,:]
    A2 = m.A2[dynamics,:]
    A3 = m.A3[dynamics,:]
    A4 = m.A4[dynamics,:]
    B1 = m.B1[dynamics,:]
    B2 = m.B2[dynamics,:]
    B3 = m.B3[dynamics,:]
    B4 = m.B4[dynamics,:]
    B5 = m.B5[dynamics,:]
    C1 = m.C1_states[dynamics,:]
    C2 = m.C2_states[dynamics,:]

    A1 += B5*staticsol.YA1
    A2 += B5*staticsol.YA2
    A3 += B5*staticsol.YA3
    A4 += B5*staticsol.YA4
    B1 += B5*staticsol.YB1
    B2 += B5*staticsol.YB2
    B3 += B5*staticsol.YB3
    B4 += B5*staticsol.YB4
    C1 += B5*staticsol.YC1
    C2 += B5*staticsol.YC2

    return rationalize.(BigInt, A1)*(model.perturbationsystem.Jx*model.perturbationsystem.Xs[:,:,r] + model.perturbationsystem.Jz*rationalize.(BigInt, exosol.Z)) + rationalize.(BigInt, A2)*(rationalize.(BigInt, Wz)*rationalize.(BigInt, exosol.Z)) + rationalize.(BigInt, A3)*(model.perturbationsystem.Xx*model.perturbationsystem.Xs[:,:,r] + model.perturbationsystem.Xz*rationalize.(BigInt, exosol.Z)) + rationalize.(BigInt, A4*(Zz*exosol.Z)) + rationalize.(BigInt, B1)*model.perturbationsystem.Js[:,:,r] + rationalize.(BigInt, B2*exosol.W) + rationalize.(BigInt, B3)*model.perturbationsystem.Xs[:,:,r] + rationalize.(BigInt, B4*exosol.Z) + rationalize.(BigInt, hcat(C1, C2))
end

function createsystem(model::RSDSGEModel,m::Vector{EvaluatedMatrices},statics::Vector{StaticSub},exos::Vector{ExoSolution},Wx::Matrix{Float64},Wz::Matrix{Float64})
    dynamics = map(((x,y),) -> !x && !y, zip(model.equations.isstatic, model.equations.isexo))
    system = Array{Symbolics.Num}(undef, count(!iszero, dynamics)*model.meta.numregimes,model.meta.numstates)
    for r in 1:model.meta.numregimes
        system[(1+(r-1)*count(!iszero, dynamics)):r*count(!iszero, dynamics),:] = createsystem(model,m[r],statics[r],exos[r],Wx,Wz,dynamics,r)
    end
    return reshape(system,length(system))
end

function markequations(model::RSDSGEModel,block)
    dynamics = .!model.equations.isstatic .& .!model.equations.isexo
    forblock = falses(0,model.meta.numstates)
    removefromsystem = falses(0,model.meta.numstates)
    keepeq = falses(count(!iszero, dynamics),model.meta.numstates)
    keepeq[block.equations[dynamics],block.states] = true
    removeeq = falses(count(!iszero, dynamics),model.meta.numstates)
    removeeq[block.equations[dynamics],:] = true
    for r in 1:model.meta.numregimes
        forblock = vcat(forblock,keepeq)
        removefromsystem = vcat(removefromsystem,removeeq)
    end
    return forblock,removefromsystem
end

function replaceblocks(model,block,system)
    # It can still depend on any or all of the exos, so we don't sub any of those out.
    tosubout = model.perturbationsystem.Xs[block.vars[model.perturbationsystem.isX],.!block.states,:]
    tosubout = reshape(tosubout,length(tosubout))
    tosubout2 = model.perturbationsystem.Js[block.vars[model.perturbationsystem.isJ],.!block.states,:]
    tosubout = vcat(tosubout,reshape(tosubout2,length(tosubout2)))

    zeroarray = Array{Symbolics.Num}(length(tosubout))
    fill!(zeroarray, Symbolics.Num(0))
    replacearg = createpairedtuplelist(tosubout,zeroarray)
    if !isempty(replacearg) > 0
        return squeeze(Symbolics.substitute(system,replacearg...),2),Dict{Symbolics.Num,Symbolics.Num}(zip(tosubout,zeroarray)) # I squeeze because subs returns a 2d with singleton trailing for some reason
    else
        return system,Dict{Symbolics.Num,Symbolics.Num}()
    end
end

function separateblocks(model::RSDSGEModel,system::Vector{Symbolics.Num})
    restricteddecisionrules = Dict{Symbolics.Num,Symbolics.Num}() # To hold the decision rule solutions that are restricted to zero because of block exogeneity
    blocks = Vector{Vector{Symbolics.Num}}()
    usedinblock = falses(length(system))
    for block in model.perturbationsystem.blocks
        blockwhich,removewhich = markequations(model,block)
        blockwhich = reshape(blockwhich,length(system))
        removewhich = reshape(removewhich,length(system))
        usedinblock = usedinblock .| removewhich
        # Now, replace any decision rules for states that don't exist in the block by zero
        blocksystem,blockdecisionrules = replaceblocks(model,block,system[blockwhich])
        push!(blocks,blocksystem)
        merge!(restricteddecisionrules,blockdecisionrules)
    end
    return blocks,system[.!usedinblock],restricteddecisionrules
end

function subsdict(target,dictionary)
    error("nope, replace me") # TODO: Remove this function
    subsargs = ()
    if !isempty(dictionary) > 0
        for (key,val) in dictionary
            subsargs = (subsargs...,(key,val))
        end
        return SymPy.subs(target,subsargs...)
    else
        return target
    end
end

function solveblocks(blocks::Vector{Vector{Symbolics.Num}},system::Vector{Symbolics.Num},restricteddecisionrules::Dict{Symbolics.Num, Symbolics.Num},verbose::Bool,method::Symbol)
    blockdict = Dict{Symbolics.Num,Float64}() # To hold the solutions to all the blocks
    if length(blocks)>0
        verbose && println("Solving the $(length(blocks)) model sub-blocks:")
        for (i,block) in enumerate(blocks)
            printif(verbose,"Block $(i) of $(length(blocks)) ($(length(block)) equations)...")
            if method == :symbolic
                sol = symbolic_solve(block)
                if length(sol) > 1
                    error("Multiple solutions for blocks not handled")
                end
                sol = first(sol)
                if hasimaginary(sol)
                    error("Imaginary solution for block, not handled")
                end
                sol = Dict{Symbolics.Num,Float64}(k => real(v) for (k,v) in sol)
            elseif method == :numeric
                sol = numericsolve(reshape(block,length(block)))
            else
                error("Option $(method) is not a valid solution method.")
            end
            merge!(blockdict,sol)
            verbose && println("done.")
        end
        verbose && print("Subbing block solutions into the full system...")
        system = subsdict(system,blockdict)
        system = subsdict(system,restricteddecisionrules)
        if ndims(system) > 1 && !isempty(system)
             # Squeeze because subs returns trailing singleton
            system = squeeze(system,2)
        end
        verbose && println("done.")
    end
    return system,blockdict
end

function convert_to_polynomial(expr::SymbolicUtils.BasicSymbolic)
    deg = Symbolics.degree(expr)
    var = Symbolics.get_variables(expr)
    if length(var) > 1
        error("More than one unknown in polynomial")
    end
    var = first(var)
    coeffs = ComplexF64[ComplexF64.(Symbolics.coeff(expr, var^d)) for d in 0:deg]
    @assert length(coeffs) == deg + 1
    return Polynomials.Polynomial(coeffs)
end

function _iterate_next_solution!(solutions::Vector{Dict{Symbolics.Num,Union{ComplexF64,Float64}}}, idx::Int, gb::Vector{<:Symbolics.BasicSymbolic}, solve_for::Symbolics.BasicSymbolic)
    this_solution = solutions[idx]
    num_solved_for = length(this_solution)
    gb_eq = gb[num_solved_for + 1]
    eq = Symbolics.substitute(gb_eq, this_solution) # Pass in what we already know
    poly = convert_to_polynomial(eq)
    roots = Polynomials.roots(poly)

    if length(roots) == 1
        # Easy, single solution, so just mutate this Solution
        solutions[idx][solve_for] = roots[1]
        return solutions
    elseif length(roots) > 1
        # Harder, need to mutate this solutoin with the first root, and add new solution
        base_sol = copy(this_solution)
        solutions[idx][solve_for] = roots[1]
        for r in roots[2:end]
            newsol = copy(base_sol)
            newsol[solve_for] = r
            push!(solutions, newsol)
        end

        return solutions
    else
        error("No roots for equation; should be impossible")
    end
end


function symbolic_solve(system::Vector{Symbolics.Num})
    if !isempty(system)
        unknowns = unique(vcat(Symbolics.get_variables.(system)...))
        if length(system) == length(unknowns)
            # Pre-allocate solution storage
            solutions = Vector{Dict{Symbolics.Num,Union{ComplexF64,Float64}}}()

            # Rewrite the system into a solveable form
            gb = Symbolics.groebner_basis(system, ordering = Groebner.Lex(unknowns...))
            # First equation will be in terms of only the first variable in unknowns (and
            # then recursively down the list - that's what the Lex ordering does)
            gb = Symbolics.simplify.(gb) # gb is in polyform from symbolicutils, which doesn't play nice
            poly = convert_to_polynomial(gb[1])# Need to use root finding from Polynomials.jl, so convert to that

            # Get the roots for the first equation to seed the solution, and then start
            # iterating (breadth-ish first) to flesh out all possible solution combinations
            roots = Polynomials.roots(poly)
            for root in roots
                push!(solutions, Dict{Symbolics.Num,Union{ComplexF64,Float64}}(unknowns[end] => root))
            end

            while any(x -> length(x) != length(unknowns), solutions)
                # We have at least one solution in the solutions array that has not been
                # fully fleshed out (ie the dictionary does not contain a solution for
                # at least one of the unknowns)
                idx = findfirst(x -> length(x) != length(unknowns), solutions)
                _iterate_next_solution!(solutions, idx, gb, unknowns[end - length(solutions[idx])])
            end

            return solutions
        else
            error("Poorly specified system. You have $(length(system)) equations for $(length(Symbolics.get_variables(system))) unknowns.")
        end
    else
        # Sub-blocks account for the whole system. This is a no-op
        return Dict{Symbolics.Num,Symbolics.Num}()
    end
end

function lambdifysystem(system::Vector{Symbolics.Num},args)
    return Function[Symbolics.build_function(expr,args...,expression=Val{false}) for expr in system]
end

function evaluatearray(x::Vector{<:Number},farray::Vector{Function})
    return map(f -> f(x...), farray)
end

function numericsolve(system::Vector{Symbolics.Num})
    if !isempty(system)
        symbols = unique(vcat(Symbolics.get_variables.(system)...))
        tosolve = lambdifysystem(system,symbols)
        fwrap = x->evaluatearray(x,tosolve)
        nlout = NLsolve.nlsolve(NLsolve.not_in_place(fwrap), zeros(length(symbols)))
        if !NLsolve.converged(nlout)
            error("Numerical solver failed to converge.")
        end
        return Dict{Symbolics.Num,Float64}(zip(symbols,nlout.zero))
    else
        # Sub-blocks account for the whole system. This is a no-op
        return Dict{Symbolics.Num,Any}()
    end
end

function hasimaginary(sol)
    for val in values(sol)
        if imag(val) > 1e-12
            return true
        end
    end
    return false
end

function formatdynamicsolution(model::RSDSGEModel,solution::Dict)
    return DynamicSolution[
        DynamicSolution(
            Symbolics.value.(map(i -> solution[i], model.perturbationsystem.Xs[:,:,r])), # have to unwrap the Symbolics.Num
            Symbolics.value.(map(i -> solution[i], model.perturbationsystem.Js[:,:,r]))
        )
        for r in 1:model.meta.numregimes
    ]
end

function arrangesolutionmatrices(model::RSDSGEModel,dynamicsolution::Vector{DynamicSolution})
    Xx = zeros(Float64, 0,count(!iszero, model.perturbationsystem.isX[model.vars.isstate]))
    Xz = zeros(Float64, 0,count(!iszero, model.perturbationsystem.isZ[model.vars.isstate]))
    Jx = zeros(Float64, 0,count(!iszero, model.perturbationsystem.isX[model.vars.isstate]))
    Jz = zeros(Float64, 0,count(!iszero, model.perturbationsystem.isZ[model.vars.isstate]))
    for r in 1:model.meta.numregimes
        Xx = vcat(Xx,dynamicsolution[r].X[:,model.perturbationsystem.isX[model.vars.isstate]])
        Jx = vcat(Jx,dynamicsolution[r].J[:,model.perturbationsystem.isX[model.vars.isstate]])
        Xz = vcat(Xz,dynamicsolution[r].X[:,model.perturbationsystem.isZ[model.vars.isstate]])
        Jz = vcat(Jz,dynamicsolution[r].J[:,model.perturbationsystem.isZ[model.vars.isstate]])
    end
    return Xx,Xz,Jx,Jz
end

function solvestatics(model::RSDSGEModel,staticsols::Vector{StaticSub},exosols::Vector{ExoSolution},dynamicsols::Vector{DynamicSolution},Xx::Matrix{Float64},Xz::Matrix{Float64},Jx::Matrix{Float64},Jz::Matrix{Float64},Wz::Matrix{Float64},Zz::Matrix{Float64})
    # This function uses the exo and dynamic solutions to return the actual solutions for static variables, rather
    # than 'solutions' that just express them as a function of other variables
    return StaticSolution[
        StaticSolution(static.YA1*(Jx*dynamic.X + Jz*exo.Z) + static.YA2*(Wz*exo.Z) + static.YA3*(Xx*dynamic.X + Xz*exo.Z) + static.YA4*(Zz*exo.Z) + static.YB1*dynamic.J + static.YB2*exo.W + static.YB3*dynamic.X + static.YB4*exo.Z + hcat(static.YC1, static.YC2))
        for (r,static,exo,dynamic) in zip(1:model.meta.numregimes,staticsols,exosols,dynamicsols)
    ]
end

function convert2decision(model,dynamics,exos,statics)
    decision = zeros(Float64,model.meta.numvars,model.meta.numvars,model.meta.numregimes)
    Xindices = model.vars.states_indices[model.perturbationsystem.isX[model.vars.isstate]]
    Zindices = model.vars.states_indices[model.perturbationsystem.isZ[model.vars.isstate]]
    indices = vcat(Xindices,Zindices)

    for r in 1:model.meta.numregimes
        decision[model.perturbationsystem.isY,indices,r] = statics[r].Y
        decision[model.perturbationsystem.isJ,indices,r] = dynamics[r].J
        decision[model.perturbationsystem.isX,indices,r] = dynamics[r].X
        decision[model.perturbationsystem.isZ,indices,r] = exos[r].Z
        decision[model.perturbationsystem.isW,indices,r] = exos[r].W
    end

    return decision
end

function solveshocks(m::EvaluatedMatrices,Xx::Matrix{Float64},Xz::Matrix{Float64},Jx::Matrix{Float64},Jz::Matrix{Float64},Wz::Matrix{Float64},Zz::Matrix{Float64})
    return ShocksSolution(-hcat(m.B1, m.B2, m.A1*Jx+m.A3*Xx+m.B3, m.A1*Jz+m.A2*Wz+m.A3*Xz+m.A4*Zz+m.B4, m.B5)\m.C_shocks)
end

function rearrangeshocks(model::RSDSGEModel,ordered::Vector{ShocksSolution})
    # J,W,X,Z,Y
    shockssol = zeros(Float64,model.meta.numvars,model.meta.numshocks,model.meta.numregimes)
    for r in 1:model.meta.numregimes
        offset = 1
        shockssol[model.perturbationsystem.isJ,:,r] = ordered[r].orderedME[offset:(offset-1+count(!iszero, model.perturbationsystem.isJ)),:]
        offset += count(!iszero, model.perturbationsystem.isJ)
        shockssol[model.perturbationsystem.isW,:,r] = ordered[r].orderedME[offset:(offset-1+count(!iszero, model.perturbationsystem.isW)),:]
        offset += count(!iszero, model.perturbationsystem.isW)
        shockssol[model.perturbationsystem.isX,:,r] = ordered[r].orderedME[offset:(offset-1+count(!iszero, model.perturbationsystem.isX)),:]
        offset += count(!iszero, model.perturbationsystem.isX)
        shockssol[model.perturbationsystem.isZ,:,r] = ordered[r].orderedME[offset:(offset-1+count(!iszero, model.perturbationsystem.isZ)),:]
        offset += count(!iszero, model.perturbationsystem.isZ)
        shockssol[model.perturbationsystem.isY,:,r] = ordered[r].orderedME[offset:(offset-1+count(!iszero, model.perturbationsystem.isY)),:]
    end
    return shockssol
end

function solveshocks(model::RSDSGEModel,m::Vector{EvaluatedMatrices},Xx::Matrix{Float64},Xz::Matrix{Float64},Jx::Matrix{Float64},Jz::Matrix{Float64},Wz::Matrix{Float64},Zz::Matrix{Float64})
    out = ShocksSolution[solveshocks(m[r],Xx,Xz,Jx,Jz,Wz,Zz) for r in 1:model.meta.numregimes]
    # Now rearrange it into the same order as the variables are listed in the model
    shockssol = rearrangeshocks(model,out)
    return shockssol
end

function solveconstant(m::EvaluatedMatrices,Xx::Matrix{Float64},Xz::Matrix{Float64},Jx::Matrix{Float64},Jz::Matrix{Float64},Wz::Matrix{Float64},Zz::Matrix{Float64})
    return ConstantSolution(-hcat(m.B1, m.B2, m.A1*Jx+m.A3*Xx+m.B3, m.A1*Jz+m.A2*Wz+m.A3*Xz+m.A4*Zz+m.B4, m.B5)\m.C_constant)
end

function rearrangeconstant(model,ordered)
    # J,W,X,Z,Y
    constantsol = zeros(Float64,model.meta.numvars,1,model.meta.numregimes)
    for r in 1:model.meta.numregimes
        offset = 1
        constantsol[model.perturbationsystem.isJ,1,r] = ordered[r].orderedMC[offset:(offset-1+count(!iszero, model.perturbationsystem.isJ)),1]
        offset += count(!iszero, model.perturbationsystem.isJ)
        constantsol[model.perturbationsystem.isW,1,r] = ordered[r].orderedMC[offset:(offset-1+count(!iszero, model.perturbationsystem.isW)),1]
        offset += count(!iszero, model.perturbationsystem.isW)
        constantsol[model.perturbationsystem.isX,1,r] = ordered[r].orderedMC[offset:(offset-1+count(!iszero, model.perturbationsystem.isX)),1]
        offset += count(!iszero, model.perturbationsystem.isX)
        constantsol[model.perturbationsystem.isZ,1,r] = ordered[r].orderedMC[offset:(offset-1+count(!iszero, model.perturbationsystem.isZ)),1]
        offset += count(!iszero, model.perturbationsystem.isZ)
        constantsol[model.perturbationsystem.isY,1,r] = ordered[r].orderedMC[offset:(offset-1+count(!iszero, model.perturbationsystem.isY)),1]
    end
    return constantsol
end

function addssvalues!(constantsol,model)
    for r in 1:model.meta.numregimes
        constantsol[:,1,r] += model.steadystate.values
    end
    return
end

function solveconstant(model::RSDSGEModel,m::Vector{EvaluatedMatrices},Xx::Matrix{Float64},Xz::Matrix{Float64},Jx::Matrix{Float64},Jz::Matrix{Float64},Wz::Matrix{Float64},Zz::Matrix{Float64})
    out = ConstantSolution[solveconstant(m[r],Xx,Xz,Jx,Jz,Wz,Zz) for r in 1:model.meta.numregimes]
    # Now rearrange it into the same order as the variables are listed in the model
    constantsol = rearrangeconstant(model,out)
    # Add steady state values
    addssvalues!(constantsol,model)
    return constantsol
end

function solve(model::RSDSGEModel,verbose::Bool = true, method::Symbol = :symbolic)
    verbose && print("Evaluating coefficient matrices...")
    m = evaluatematrices(model)
    verbose && println("done.")

    verbose && print("Subbing out static variables...")
    statics = substatics(model,m)
    verbose && println("done.")

    verbose && print("Solving exogenous equations...")
    exosolutions = solveexos(model,m,statics)
    Wz,Zz = arrangeexomatrices(model,exosolutions)
    verbose && println("done.")

    verbose && print("Constructing system of quadratic equations...")
    system = createsystem(model,m,statics,exosolutions,Wz,Zz)
    verbose && println("done")

    blocks,system,restricteddecisionrules = separateblocks(model,system)
    system,blockdict = solveblocks(blocks,system,restricteddecisionrules,verbose,method)
    merge!(blockdict,restricteddecisionrules)

    verbose && print("Solving complete system of quadratic equations ($(length(system)) equations)...")
    dynamicsolutions = if method == :symbolic
        symbolic_solve(system)
    elseif methods == :numeric
        numericsolve(system)
    else
        error("Unknown solution method $method")
    end

    verbose && println("done")

    if isa(dynamicsolutions,Dict)
        dynamicsolutions = [dynamicsolutions]
    end
    println("Found $(length(dynamicsolutions)) solution(s).")
    solutions = Vector{RSDSGESolution}()
    for (i,sol) in enumerate(dynamicsolutions)
        # Ignore any imaginary solutions
        if !hasimaginary(sol) # This tests approximate zero for imaginary part
            # Strip out the imaginary parts
            # Merge in the block dict solutions
            sol = Dict{Symbolics.Num,Float64}(k => real(v) for (k,v) in sol)
            sol = merge(sol,blockdict)
            # TODO: Here's where I probably want to test for stability using FRWZ MSS method?

            # Format into nicer form (i.e., not a dictionary)
            dynamics = formatdynamicsolution(model,sol)
            # Arrange the states solution in Xx, Xz etc
            Xx,Xz,Jx,Jz = arrangesolutionmatrices(model,dynamics)
            # Back out the statics
            staticsol = solvestatics(model,statics,exosolutions,dynamics,Xx,Xz,Jx,Jz,Wz,Zz)
            # Format into decision rule matrix
            decisions = convert2decision(model,dynamics,exosolutions,staticsol)
            # Solve shocks
            shocks = solveshocks(model,m,Xx,Xz,Jx,Jz,Wz,Zz)
            # Solve constant
            constant = solveconstant(model,m,Xx,Xz,Jx,Jz,Wz,Zz)
            # Wrap into solution type and push into array
            push!(solutions,RSDSGESolution(decisions,shocks,constant))
        else
            println("Solution number $(i) has imaginary components, ignoring...")
        end
    end
    return solutions
end
