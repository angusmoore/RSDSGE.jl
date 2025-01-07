function createarglist(parameters,vars,transmatrix)
    arg = Array{Symbolics.Num}(undef, 1)
    arg = [transmatrix.completelygeneric] # Transition prob first
    append!(arg,vars.contemps_sym) # Add in the steady state values
    # Non-switching vars first
    append!(arg,parameters.generic_sym[1,.!parameters.isswitching])
    # For parameters that affect SS, just need the bar version
    append!(arg,parameters.bar_sym[parameters.affectsSS])
    # Finally, for parameters that switch but _don't_ affect SS, need both R and and RP version
    append!(arg,parameters.generic_sym[1,parameters.isswitching .& .!parameters.affectsSS])
    append!(arg,parameters.generic_sym[2,parameters.isswitching .& .!parameters.affectsSS])
    return arg
end

function lambdifycoefficient(parameters,vars,transmatrix,symmatrix,arglist)
    fmatrix = Array{Function}(undef,size(symmatrix)...)
    for i in eachindex(symmatrix)
        if symmatrix[i] isa Symbolics.Num
            fmatrix[i] = Symbolics.build_function(symmatrix[i],arglist...,expression=Val{false})
        else
            fmatrix[i] = x -> 0
        end
    end
    return fmatrix
end

function createXx(meta,vars,Xs,usewhichstates)
    select = usewhichstates[vars.isstate]
    allregimes = Xs[:,select,:]
    out = Array{Symbolics.Num}(undef,0,count(!iszero, usewhichstates))
    for r = 1:meta.numregimes
        out = vcat(out,allregimes[:,:,r])
    end
    return out
end

function createXs(meta,vars,usevar)
    out = Array{Symbolics.Num}(undef,count(!iszero, usevar),meta.numstates,meta.numregimes)
    for r = 1:meta.numregimes
        for (sstep,sindex) in enumerate(vars.states_indices)
            statename = vars.names[sindex]
            for (i,varname) in enumerate(vars.names[usevar])
                out[i,sstep,r] = Symbolics.variable("$(varname)_diff_$(statename)_R$(r)")
            end
        end
    end
    return out
end

function createcoefficient(equations::Array{Symbolics.Num, 1},vars::Array{Symbolics.Num, 1},transprob,subsargs::AbstractDict)
    deriv = Array{Symbolics.Num}(undef, length(equations), length(vars))
    for (i, equation) in enumerate(equations)
        for (j, var) in enumerate(vars)
            deriv[i, j] = Symbolics.derivative(equation, var)
        end
    end
    deriv = Symbolics.substitute(deriv,subsargs)
    return transprob*deriv
end

function createleadcoefficients(subsargs::AbstractDict,transmatrix,equations,jump,nonstateexo,endo,exo)
    # These are all for one specific future regime, so I use generic transition and parameters
    A1 = createcoefficient(equations,jump,transmatrix.completelygeneric,subsargs)
    A2 = createcoefficient(equations,nonstateexo,transmatrix.completelygeneric,subsargs)
    A3 = createcoefficient(equations,endo,transmatrix.completelygeneric,subsargs)
    A4 = createcoefficient(equations,exo,transmatrix.completelygeneric,subsargs)

    return A1,A2,A3,A4
end

function createcontempcoefficients(subsargs,transmatrix,equations,jump,nonstateexo,endo,exo,static)
    # These are all for one specific future regime, so I use generic transition and parameters
    B1 = createcoefficient(equations,jump,transmatrix.completelygeneric,subsargs)
    B2 = createcoefficient(equations,nonstateexo,transmatrix.completelygeneric,subsargs)
    B3 = createcoefficient(equations,endo,transmatrix.completelygeneric,subsargs)
    B4 = createcoefficient(equations,exo,transmatrix.completelygeneric,subsargs)
    B5 = createcoefficient(equations,static,transmatrix.completelygeneric,subsargs)

    return B1,B2,B3,B4,B5
end

function createstatescoefficients(subsargs,transmatrix,equations,endo,exo)
    C1 = createcoefficient(equations,endo,transmatrix.completelygeneric,subsargs)
    C2 = createcoefficient(equations,exo,transmatrix.completelygeneric,subsargs)

    return C1,C2
end

function createshockscoefficients(subsargs,transmatrix,equations,shocks)
    C_shocks = createcoefficient(equations,shocks,transmatrix.completelygeneric,subsargs)
    return C_shocks
end

function createparametercoefficients(subsargs,transmatrix,equations,parameters)
    C_constant = createcoefficient(equations,parameters,transmatrix.completelygeneric,subsargs)
    return C_constant
end

function lambdifyall(A1,A2,A3,A4,B1,B2,B3,B4,B5,C1,C2,C_shocks,C_constant,C_constant_p,parameters,vars,transmatrix)
    # Construct the arg list first
    arglist = createarglist(parameters,vars,transmatrix)
    A1_L = lambdifycoefficient(parameters,vars,transmatrix,A1,arglist)
    A2_L = lambdifycoefficient(parameters,vars,transmatrix,A2,arglist)
    A3_L = lambdifycoefficient(parameters,vars,transmatrix,A3,arglist)
    A4_L = lambdifycoefficient(parameters,vars,transmatrix,A4,arglist)
    B1_L = lambdifycoefficient(parameters,vars,transmatrix,B1,arglist)
    B2_L = lambdifycoefficient(parameters,vars,transmatrix,B2,arglist)
    B3_L = lambdifycoefficient(parameters,vars,transmatrix,B3,arglist)
    B4_L = lambdifycoefficient(parameters,vars,transmatrix,B4,arglist)
    B5_L = lambdifycoefficient(parameters,vars,transmatrix,B5,arglist)
    C1_L = lambdifycoefficient(parameters,vars,transmatrix,C1,arglist)
    C2_L = lambdifycoefficient(parameters,vars,transmatrix,C2,arglist)
    C_shocks_L = lambdifycoefficient(parameters,vars,transmatrix,C_shocks,arglist)
    C_constant_L = lambdifycoefficient(parameters,vars,transmatrix,C_constant,arglist)
    C_constant_p_L = lambdifycoefficient(parameters,vars,transmatrix,C_constant_p,arglist)

    return A1_L,A2_L,A3_L,A4_L,B1_L,B2_L,B3_L,B4_L,B5_L,C1_L,C2_L,C_shocks_L,C_constant_L,C_constant_p_L
end

function comboasbitarray(i,numeqs)
    return digits(Bool, i, base = 2, pad = ceil(Int,numeqs))
end

function whichvariables(vars,equation)
    # Find out which vars are in the equations
    has_leads = expressioncontains(equation,vars.leads_sym)
    has_contemps = expressioncontains(equation,vars.contemps_sym)
    has_state = expressioncontains(equation, vars.states_sym)
    # now union across and negate out any exogenous vars
    has_vars = (has_leads .| has_contemps)
    # Make any states true as well
    has_vars[vars.states_indices[has_state]] .= true

    # Union across any contemp or lead versions of the state
    has_state = has_vars[vars.states_indices] .| has_state

    return has_vars, has_state
end

function getblocksize(block::Block)
    return block.size
end

function checkredundancy(vars,block,list)
    possibilities = length(list)
    numcombos = 2^possibilities-1
    # Do the easiest check first. Is it just an existing block with tacked on exo or static
    for j in 1:length(list)
        combo = list[j].vars
        if all(vars.isexo[block.vars. & .!combo] .| vars.isstatic[block.vars .& .!combo])
            return true
        end
    end

    for i = 1:numcombos
        useblocks = comboasbitarray(i,possibilities)
        combo = falses(length(list[1].vars))
        for j in 1:length(useblocks)
            if useblocks[j]
                combo = combo .| list[j].vars
            end
        end
        if combo == block.vars
            return true
        end
        # Finally, if all the block is doing is adding an extra exogenous equation, that's not useful either
        if all(vars.isexo[block.vars .& .!combo])
            return true
        end
    end
    return false
end

function removeredundant(vars,blocks)
    isredundant = falses(length(blocks))
    for i = length(blocks):-1:1
        isredundant[i] = checkredundancy(vars,blocks[i],blocks[1:(i-1)])
    end
    println(" $(count(!iszero, isredundant)) of them are redundant.")
    return blocks[.!isredundant]
end

function findblocks(meta,vars,equations)
    # First, figure out which variables are in which equations
    numeqs = meta.numvars
    expressionsymbols = fill(false, meta.numvars,numeqs)
    expressionstates = fill(false, meta.numstates,numeqs)
    for (i,equation) in enumerate(equations.syms)
        expressionsymbols[:,i],expressionstates[:,i] = whichvariables(vars,equation)
    end

    numcombos = 2^numeqs-1
    blocks = Array{Block}(undef, 0)
    for i = 1:numcombos
        whichequations = comboasbitarray(i,numeqs)
        blocknumeqs = count(!iszero, whichequations)
        blockvars = any(expressionsymbols[:,whichequations], dims = 2)
        blockvars = reshape(blockvars, length(blockvars))
        blocknumvars = count(!iszero, blockvars)
        if blocknumeqs == blocknumvars && blocknumeqs != numeqs
            # Candidate blocks. First check if it's all exogenous. If so, don't bother.
            # Similarly, if the block is just all equations except statics or exos, don't bother.
            # Finally, if it contains more than 60 per cent of the equations, don't bother. This is an ad hoc way to avoid a very large number of blocks, which enourmously slows down the code.
            if !all(vars.isexo[blockvars]) && !all(equations.isstatic[.!whichequations] .| equations.isexo[.!whichequations]) && blocknumeqs < 0.6*numeqs
                push!(blocks,Block(blocknumvars,whichequations,blockvars,squeeze(any(expressionstates[:,whichequations],2),2)))
            end
        end
    end
    sort!(blocks,by=getblocksize)
    print("Found $(length(blocks)) candidate blocks...")
    return removeredundant(vars,blocks)
end

function PerturbationSystem(meta,parameters,vars,shocks,equations,steadystate,transmatrix)

    """
    Variable classifications:
    J: Jump
    X: Endogenous state
    W: Exogenous non-state
    Z: Exogenous state
    Y: static

    A1: d/dJ(+1)
    A2: d/dW(+1)
    A3: d/dX(+1)
    A4: d/dZ(+1)

    B1: d/dJ
    B2: d/dW
    B3: d/dX
    B4: d/dZ
    B5: d/dY

    C3_state: d/dX(-1)
    C4_state: d/dZ(-1)
    C_shock: d/de
    C_constant: d/dPERTURBATION

    """

    # First, figure out the classification of each variable
    isJ = .!vars.isstate .& .!vars.isexo .& .!vars.isstatic
    isY = copy(vars.isstatic)
    isW = .!vars.isstate .& vars.isexo
    isZ = vars.isstate .& vars.isexo
    isX = vars.isstate .& .!vars.isexo
    varnums = range(1,meta.numvars)
    Jindices = varnums[isJ]
    Yindices = varnums[isY]
    Windices = varnums[isW]
    Zindices = varnums[isZ]
    Xindices = varnums[isX]

    println("----- Diagnostics -----")
    println("Jump variables: $(count(!iszero, isJ))")
    println("Static variables: $(count(!iszero, isY))")
    println("Non-state exo variables: $(count(!iszero, isW))")
    println("Exo states variables: $(count(!iszero, isZ))")
    println("Endo states variables: $(count(!iszero, isX))")

    # Create the matrices of Symbolic taylor unknowns, which are used by solve. We don't need them for exos or statics, since
    # we solve for these linearly/sub them out (respectively)

    # Generic function to create a 3-d array of specific vars (3rd arg) wrt all state variables.
    # the third dimension of the array indexes regimes
    Xs = createXs(meta,vars,isX)
    Js = createXs(meta,vars,isJ)

    # Generic function to stack Xs into a column of all states. Restricts which states to keep the derivatives of (3rd argument)
    Xx = createXx(meta,vars,Xs,isX)
    Xz = createXx(meta,vars,Xs,isZ)
    Jx = createXx(meta,vars,Js,isX)
    Jz = createXx(meta,vars,Js,isZ)

    subsargs = Dict(
        old => new
        for (old, new) in zip(
            [vars.leads_sym; vars.states_sym; shocks.syms; parameters.generic_sym[1,parameters.affectsSS]; parameters.generic_sym[2,parameters.affectsSS]],
            [vars.contemps_sym; vars.contemps_sym[vars.states_indices]; zeros(Float64,length(shocks.syms)); parameters.bar_sym[parameters.affectsSS]; parameters.bar_sym[parameters.affectsSS]; parameters.bar_sym[parameters.affectsSS]]
        )
    )

    A1,A2,A3,A4 = createleadcoefficients(subsargs,transmatrix,equations.syms,vars.leads_sym[isJ],vars.leads_sym[isW],vars.leads_sym[isX],vars.leads_sym[isZ])
    B1,B2,B3,B4,B5 = createcontempcoefficients(subsargs,transmatrix,equations.syms,vars.contemps_sym[isJ],vars.contemps_sym[isW],vars.contemps_sym[isX],vars.contemps_sym[isZ],vars.contemps_sym[isY])
    C1,C2 = createstatescoefficients(subsargs,transmatrix,equations.syms,vars.states_sym[isX[vars.isstate]],vars.states_sym[isZ[vars.isstate]])
    C_shocks = createshockscoefficients(subsargs,transmatrix,equations.syms,shocks.syms)
    C_constant = createparametercoefficients(subsargs,transmatrix,equations.syms,parameters.generic_sym[1,parameters.affectsSS])
    C_constant_p = createparametercoefficients(subsargs,transmatrix,equations.syms,parameters.generic_sym[2,parameters.affectsSS])

    A1_L,A2_L,A3_L,A4_L,B1_L,B2_L,B3_L,B4_L,B5_L,C1_L,C2_L,C_shocks_L,C_constant_L,C_constant_p_L = lambdifyall(A1,A2,A3,A4,B1,B2,B3,B4,B5,C1,C2,C_shocks,C_constant,C_constant_p,parameters,vars,transmatrix)

    blocks = findblocks(meta,vars,equations)

    return PerturbationSystem(A1_L,A2_L,A3_L,A4_L,B1_L,B2_L,B3_L,B4_L,B5_L,C1_L,C2_L,C_shocks_L,C_constant_L,C_constant_p_L,Xs,Xx,Xz,Js,Jx,Jz,isJ,isY,isW,isZ,isX,blocks)
end
