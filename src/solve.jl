type EvaluatedMatrices
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

type StaticSub
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

type ExoSolution
    W::Array{Float64,2}
    Z::Array{Float64,2}
end

type DynamicSolution
    X::Array{Float64,2}
    J::Array{Float64,2}
end

type StaticSolution
    Y::Array{Float64,2}
end

type ShocksSolution
    orderedME::Array{Float64,2}
end

type ConstantSolution
    orderedMC::Array{Float64,2}
end

type RSDSGESolution
    MX::Array{Float64,3}
    ME::Array{Float64,3}
    MC::Array{Float64,3}
end

function EvaluatedMatrices(model,r)
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
    C_constant = sum(C_constant,2)
    return EvaluatedMatrices(A1,A2,A3,A4,B1,B2,B3,B4,B5,C1,C2,C_shocks,C_constant)
end

function evaluatematrices(model)
    out = Array(EvaluatedMatrices,0)
    for r in 1:model.meta.numregimes
        push!(out,EvaluatedMatrices(model,r))
    end
    return out
end

function createarglistvalues(model,r,rp)
    arg = Array(Float64,1)
    arg = [model.transmatrix.values[r,rp]] # Transition prob first
    append!(arg,model.steadystate.values) # Add in the steady state values
    # Non-switching vars first
    append!(arg,model.parameters.values_withbar[r,!model.parameters.isswitching])
    # For parameters that affect SS, just need the bar version
    append!(arg,model.parameters.values_withbar[r,model.parameters.affectsSS])
    # Finally, for parameters that switch but _don't_ affect SS, need both R and and RP version
    append!(arg,model.parameters.values_withbar[r,model.parameters.isswitching & !model.parameters.affectsSS])
    append!(arg,model.parameters.values_withbar[rp,model.parameters.isswitching & !model.parameters.affectsSS])
    return arg
end

function evaluate(model,fmatrix,r,rp)
    args = createarglistvalues(model,r,rp)
    out = zeros(Float64,size(fmatrix)...)
    if nnz(fmatrix) > 0
        R,C,F = findnz(fmatrix)
        for (r,c,f) in zip(R,C,F)
            out[r,c] = f(args...)
        end
    end
    return out
end

function createsummed(model,fmatrix,r)
    out = zeros(Float64,size(fmatrix)...)
    if !isempty(fmatrix)
        for rp in 1:model.meta.numregimes
            out += evaluate(model,fmatrix,r,rp)
        end
    end
    return out
end

function createstacked(model,fmatrix,r)
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
                x = (rp-1)*countnz(model.parameters.affectsSS) + i
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

function substatics(model,m,r)
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

function substatics(model,m)
    solutions = Array(StaticSub,model.meta.numregimes)
    for r in 1:model.meta.numregimes
        solutions[r] = substatics(model,m[r],r)
    end
    return solutions
end

function solveexos(model,m,staticsol,r)
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
    
    # These should all be non zero. This is just double checking that. I'll remove these ops for production, since
    # they carry performance cost.
    if countnz(A1)>0
        error("Non-zero A1 entries in exo block")
    end
    if countnz(A2)>0
        error("Non-zero A2 entries in exo block")
    end
    if countnz(A3)>0
        error("Non-zero A3 entries in exo block")
    end
    if countnz(A4)>0
        error("Non-zero A4 entries in exo block")
    end
    if countnz(B1)>0
        error("Non-zero B1 entries in exo block")
    end
    if countnz(B3)>0
        println(full(B3))
        error("Non-zero B3 entries in exo block")
    end

    sol = -hcat(B2,B4)\hcat(C1,C2)
    W = sol[1:countnz(model.perturbationsystem.isW),:]
    Z = sol[(1+countnz(model.perturbationsystem.isW)):end,:]
    return ExoSolution(W,Z)
end

function solveexos(model,m,staticsubs)
    solutions = Array(ExoSolution,model.meta.numregimes)
    for r in 1:model.meta.numregimes
        solutions[r] = solveexos(model,m[r],staticsubs[r],r)
    end
    return solutions
end

function arrangeexomatrices(model,exosols)
    nz = countnz(model.perturbationsystem.isZ)
    nw = countnz(model.perturbationsystem.isW)
    nx = countnz(model.perturbationsystem.isX)
    Wz = zeros(Float64,0,nz)
    Zz = zeros(Float64,0,nz)
    for r in 1:model.meta.numregimes
        Wz = vcat(Wz,exosols[r].W[:,(1+nx):end])
        Zz = vcat(Zz,exosols[r].Z[:,(1+nx):end])
    end
    return Wz,Zz
end

function createsystem(model,m,staticsol,exosol,Wz,Zz,dynamics,r)
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

    return  A1*(model.perturbationsystem.Jx*model.perturbationsystem.Xs[:,:,r] + model.perturbationsystem.Jz*exosol.Z) + A2*(Wz*exosol.Z) + A3*(model.perturbationsystem.Xx*model.perturbationsystem.Xs[:,:,r] + model.perturbationsystem.Xz*exosol.Z) + A4*(Zz*exosol.Z) + B1*model.perturbationsystem.Js[:,:,r] + B2*exosol.W + B3*model.perturbationsystem.Xs[:,:,r] + B4*exosol.Z + hcat(C1, C2) 
end

function createsystem(model,m,statics,exos,Wx,Wz)
    dynamics = !model.equations.isstatic & !model.equations.isexo
    system = Array(SymPy.Sym,countnz(dynamics)*model.meta.numregimes,model.meta.numstates)
    for r in 1:model.meta.numregimes
        system[(1+(r-1)*countnz(dynamics)):r*countnz(dynamics),:] = createsystem(model,m[r],statics[r],exos[r],Wx,Wz,dynamics,r)
    end
    return reshape(system,length(system)) # The reshape is needed because sympy can only take 1d arrays. It doesn't change anything important.
end

function markequations(model,block)
    dynamics = !model.equations.isstatic & !model.equations.isexo
    forblock = falses(0,model.meta.numstates)
    removefromsystem = falses(0,model.meta.numstates)
    keepeq = falses(countnz(dynamics),model.meta.numstates)
    keepeq[block.equations[dynamics],block.states] = true
    removeeq = falses(countnz(dynamics),model.meta.numstates)
    removeeq[block.equations[dynamics],:] = true
    for r in 1:model.meta.numregimes
        forblock = vcat(forblock,keepeq)
        removefromsystem = vcat(removefromsystem,removeeq)
    end
    return forblock,removefromsystem
end

function replaceblocks(model,block,system)
    # It can still depend on any or all of the exos, so we don't sub any of those out.
    tosubout = model.perturbationsystem.Xs[block.vars[model.perturbationsystem.isX],!block.states,:]
    tosubout = reshape(tosubout,length(tosubout))
    tosubout2 = model.perturbationsystem.Js[block.vars[model.perturbationsystem.isJ],!block.states,:]
    tosubout = vcat(tosubout,reshape(tosubout2,length(tosubout2)))
    
    sympyzero = Sym(0.0)
    zeroarray = Array(SymPy.Sym,length(tosubout))
    fill!(zeroarray,sympyzero)
    replacearg = createpairedtuplelist(tosubout,zeroarray)
    if !isempty(replacearg)>0
        return squeeze(SymPy.subs(system,replacearg...),2),Dict{SymPy.Sym,SymPy.Sym}(zip(tosubout,zeroarray)) # I squeeze because subs returns a 2d with singleton trailing for some reason
    else
        return system,Dict{SymPy.Sym,SymPy.Sym}()
    end
end

function separateblocks(model,system)
    restricteddecisionrules = Dict{SymPy.Sym,SymPy.Sym}() # To hold the decision rule solutions that are restricted to zero because of block exogeneity
    blocks = Array(Array{SymPy.Sym,1},0)
    usedinblock = falses(length(system))
    for block in model.perturbationsystem.blocks
        blockwhich,removewhich = markequations(model,block)
        blockwhich = reshape(blockwhich,length(system))
        removewhich = reshape(removewhich,length(system))
        usedinblock = usedinblock | removewhich
        # Now, replace any decision rules for states that don't exist in the block by zero
        blocksystem,blockdecisionrules = replaceblocks(model,block,system[blockwhich])
        push!(blocks,blocksystem)
        merge!(restricteddecisionrules,blockdecisionrules)
    end
    return blocks,system[!usedinblock],restricteddecisionrules
end

function subsdict(target,dictionary)
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

function solveblocks(model,blocks,system,restricteddecisionrules,verbose::Bool,method)
    blockdict = Dict{SymPy.Sym,SymPy.Sym}() # To hold the solutions to all the blocks
    if length(blocks)>0
        printlnif(verbose,"Solving the $(length(blocks)) model sub-blocks:")
        for (i,block) in enumerate(blocks)
            printif(verbose,"Block $(i) of $(length(blocks)) ($(length(block)) equations)...")
            if method == 1
                sol = solve_withfallback(block)
            elseif method == 2
                sol = numericsolve(reshape(block,length(block)))
            else
                error("Option $(method) is not a valid solution method.")
            end
            if !isa(sol,Dict)
                error("Handling multiple solutions for blocks not yet implemented.")
            end
            merge!(blockdict,sol)
            printlnif(verbose,"done.")
        end
        printif(verbose,"Subbing block solutions into the full system...")
        system = subsdict(system,blockdict)
        system = subsdict(system,restricteddecisionrules)
        if ndims(system) > 1 && !isempty(system)
             # Squeeze because subs returns trailing singleton
            system = squeeze(system,2)
        end
        printlnif(verbose,"done.")
    end
    return system,blockdict
end

function convert2rationals(system)
    for (i,expression) in enumerate(system)
        system[i]=SymPy.nsimplify(expression)
    end
    return system
end

function systemissolvable(system)
    if length(system) == length(SymPy.free_symbols(system))
        return true
    else
        println(system)
        println(SymPy.free_symbols(system))
        error("Poorly specified system. You have $(length(system)) equations for $(length(SymPy.free_symbols(system))) unknowns.")
    end
end

function solve_withfallback(system)
    # This function wraps solve to first try solving using floats. If that fails because sympy is sometimes inaccurate with floats
    # it then converts to rationals and tries again.
    if !isempty(system)
        if systemissolvable(system)
            try 
                solution = SymPy.solve(system,method="f5b")
                return solution
            catch
                warn("Solve failed. Converting to rationals and trying again.") 
                system = convert2rationals(system)
                solution = SymPy.solve(system,method="f5b")
                return solution
            end
        end
    else
        # Sub-blocks account for the whole system. This is a no-op
        return Dict{SymPy.Sym,SymPy.Sym}()
    end
end

function lambdifysystem(system,args)
    out = Array(Function,length(system))
    for (i,expr) in enumerate(system)
        out[i] = SymPy.lambdify(expr,args)
    end
    return out
end

function evaluatearray(x,farray)
    out = similar(x)
    for (i,f) in enumerate(farray)
        out[i] = f(x...)
    end
    return out
end

function numericsolve(system)
    if !isempty(system)
        symbols = SymPy.free_symbols(system)
        tosolve = lambdifysystem(system,symbols)
        fwrap = x->evaluatearray(x,tosolve)
        nlout = NLsolve.nlsolve(NLsolve.not_in_place(fwrap), zeros(length(symbols)))
        if !NLsolve.converged(nlout)
            error("Numerical solver failed to converge.")
        end
        dynamicsolutions = Dict{SymPy.Sym,Float64}(zip(symbols,nlout.zero))
        return dynamicsolutions
    else
        # Sub-blocks account for the whole system. This is a no-op
        return Dict{SymPy.Sym,SymPy.Sym}()
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

function dict2array(dict,indices)
    out = Array(valtype(dict),size(indices))
    for i in eachindex(indices)
        out[i] = dict[indices[i]]
    end
    return out    
end

function formatdynamicsolution(model,solution::Dict)
    out = Array(DynamicSolution,model.meta.numregimes)
    for r in 1:model.meta.numregimes
        Xsols = dict2array(solution,model.perturbationsystem.Xs[:,:,r])
        Jsols = dict2array(solution,model.perturbationsystem.Js[:,:,r])
        out[r] = DynamicSolution(Xsols,Jsols)
    end
    return out
end

function arrangesolutionmatrices(model,exosolution,dynamicsolution)
    Xx = Array(Float64,0,countnz(model.perturbationsystem.isX[model.vars.isstate]))
    Xz = Array(Float64,0,countnz(model.perturbationsystem.isZ[model.vars.isstate]))
    Jx = Array(Float64,0,countnz(model.perturbationsystem.isX[model.vars.isstate]))
    Jz = Array(Float64,0,countnz(model.perturbationsystem.isZ[model.vars.isstate]))
    for r in 1:model.meta.numregimes
        Xx = vcat(Xx,dynamicsolution[r].X[:,model.perturbationsystem.isX[model.vars.isstate]])
        Jx = vcat(Jx,dynamicsolution[r].J[:,model.perturbationsystem.isX[model.vars.isstate]])
        Xz = vcat(Xz,dynamicsolution[r].X[:,model.perturbationsystem.isZ[model.vars.isstate]])
        Jz = vcat(Jz,dynamicsolution[r].J[:,model.perturbationsystem.isZ[model.vars.isstate]])
    end
    return Xx,Xz,Jx,Jz
end

function solvestatics(model,staticsols,exosols,dynamicsols,Xx,Xz,Jx,Jz,Wz,Zz)
    # This function uses the exo and dynamic solutions to return the actual solutions for static variables, rather
    # than 'solutions' that just express them as a function of other variables
    out = Array(StaticSolution,model.meta.numregimes)
    for (r,static,exo,dynamic) in zip(1:model.meta.numregimes,staticsols,exosols,dynamicsols)
        out[r] = StaticSolution(static.YA1*(Jx*dynamic.X + Jz*exo.Z) + static.YA2*(Wz*exo.Z) + static.YA3*(Xx*dynamic.X + Xz*exo.Z) + static.YA4*(Zz*exo.Z) + static.YB1*dynamic.J + static.YB2*exo.W + static.YB3*dynamic.X + static.YB4*exo.Z + hcat(static.YC1, static.YC2))
    end
    return out
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

function solveshocks(model,m,Xx,Xz,Jx,Jz,Wz,Zz,r)
    return ShocksSolution(-hcat(m.B1, m.B2, m.A1*Jx+m.A3*Xx+m.B3, m.A1*Jz+m.A2*Wz+m.A3*Xz+m.A4*Zz+m.B4, m.B5)\m.C_shocks)
end

function rearrangeshocks(model,ordered)
    # J,W,X,Z,Y
    shockssol = zeros(Float64,model.meta.numvars,model.meta.numshocks,model.meta.numregimes)
    for r in 1:model.meta.numregimes
        offset = 1
        shockssol[model.perturbationsystem.isJ,:,r] = ordered[r].orderedME[offset:(offset-1+countnz(model.perturbationsystem.isJ)),:]
        offset += countnz(model.perturbationsystem.isJ)
        shockssol[model.perturbationsystem.isW,:,r] = ordered[r].orderedME[offset:(offset-1+countnz(model.perturbationsystem.isW)),:]
        offset += countnz(model.perturbationsystem.isW)
        shockssol[model.perturbationsystem.isX,:,r] = ordered[r].orderedME[offset:(offset-1+countnz(model.perturbationsystem.isX)),:]
        offset += countnz(model.perturbationsystem.isX)
        shockssol[model.perturbationsystem.isZ,:,r] = ordered[r].orderedME[offset:(offset-1+countnz(model.perturbationsystem.isZ)),:]
        offset += countnz(model.perturbationsystem.isZ)
        shockssol[model.perturbationsystem.isY,:,r] = ordered[r].orderedME[offset:(offset-1+countnz(model.perturbationsystem.isY)),:]
    end
    return shockssol
end

function solveshocks(model,m,Xx,Xz,Jx,Jz,Wz,Zz)
    out = Array(ShocksSolution,0)
    for r in 1:model.meta.numregimes
        push!(out,solveshocks(model,m[r],Xx,Xz,Jx,Jz,Wz,Zz,r))
    end
    # Now rearrange it into the same order as the variables are listed in the model
    shockssol = rearrangeshocks(model,out)
    return shockssol
end

function solveconstant(model,m,Xx,Xz,Jx,Jz,Wz,Zz,r)
    return ConstantSolution(-hcat(m.B1, m.B2, m.A1*Jx+m.A3*Xx+m.B3, m.A1*Jz+m.A2*Wz+m.A3*Xz+m.A4*Zz+m.B4, m.B5)\m.C_constant)
end

function rearrangeconstant(model,ordered)
    # J,W,X,Z,Y
    constantsol = zeros(Float64,model.meta.numvars,1,model.meta.numregimes)
    for r in 1:model.meta.numregimes
        offset = 1
        constantsol[model.perturbationsystem.isJ,1,r] = ordered[r].orderedMC[offset:(offset-1+countnz(model.perturbationsystem.isJ)),1]
        offset += countnz(model.perturbationsystem.isJ)
        constantsol[model.perturbationsystem.isW,1,r] = ordered[r].orderedMC[offset:(offset-1+countnz(model.perturbationsystem.isW)),1]
        offset += countnz(model.perturbationsystem.isW)
        constantsol[model.perturbationsystem.isX,1,r] = ordered[r].orderedMC[offset:(offset-1+countnz(model.perturbationsystem.isX)),1]
        offset += countnz(model.perturbationsystem.isX)
        constantsol[model.perturbationsystem.isZ,1,r] = ordered[r].orderedMC[offset:(offset-1+countnz(model.perturbationsystem.isZ)),1]
        offset += countnz(model.perturbationsystem.isZ)
        constantsol[model.perturbationsystem.isY,1,r] = ordered[r].orderedMC[offset:(offset-1+countnz(model.perturbationsystem.isY)),1]
    end
    return constantsol
end

function addssvalues!(constantsol,model)
    for r in 1:model.meta.numregimes
        constantsol[:,1,r] += model.steadystate.values
    end
    return
end

function solveconstant(model,m,Xx,Xz,Jx,Jz,Wz,Zz)
    out = Array(ConstantSolution,0)
    for r in 1:model.meta.numregimes
        push!(out,solveconstant(model,m[r],Xx,Xz,Jx,Jz,Wz,Zz,r))
    end
    # Now rearrange it into the same order as the variables are listed in the model
    constantsol = rearrangeconstant(model,out)
    # Add steady state values
    addssvalues!(constantsol,model)
    return constantsol
end

function printif(verbose,string)
    if verbose
        print(string)
    end
    return
end

function printlnif(verbose,string)
    if verbose
        println(string)
    end
    return
end

function solve(model::RSDSGEModel,verbose::Bool = true; method::Int=1)
    printif(verbose,"Evaluating coefficient matrices...")
    m = evaluatematrices(model)
    printlnif(verbose,"done.")
    
    printif(verbose,"Subbing out static variables...")
    statics = substatics(model,m)
    printlnif(verbose,"done.")
    
    printif(verbose,"Solving exogenous equations...")
    exosolutions = solveexos(model,m,statics)
    Wz,Zz = arrangeexomatrices(model,exosolutions)
    printlnif(verbose,"done.")
    
    printif(verbose,"Constructing system of quadratic equations...")
    system = createsystem(model,m,statics,exosolutions,Wz,Zz)
    printlnif(verbose,"done")
    
    blocks,system,restricteddecisionrules = separateblocks(model,system)
    system,blockdict = solveblocks(model,blocks,system,restricteddecisionrules,verbose,method)
    merge!(blockdict,restricteddecisionrules)

    printif(verbose,"Solving complete system of quadratic equations ($(length(system)) equations)...")
    if method == 1
        dynamicsolutions = solve_withfallback(system)
    elseif method == 2
        dynamicsolutions = numericsolve(system)
    else
        error("Option $(method) is not a valid solution method.")
    end
        
    printlnif(verbose,"done")
    
    if isa(dynamicsolutions,Dict)
        dynamicsolutions = [dynamicsolutions]
    end
    println("Found $(length(dynamicsolutions)) solution(s).")
    solutions = Array(RSDSGESolution,0)
    for (i,sol) in enumerate(dynamicsolutions)
        # Ignore any imaginary solutions
        if !hasimaginary(sol)
            # Merge in the block dict solutions
            sol = merge(sol,blockdict)
            # Here's where I probably want to test for stability using FRWZ MSS method?
            
            # Format into nicer form (i.e., not a dictionary)
            dynamics = formatdynamicsolution(model,sol)
            # Arrange the states solution in Xx, Xz etc
            Xx,Xz,Jx,Jz = arrangesolutionmatrices(model,exosolutions,dynamics) 
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
