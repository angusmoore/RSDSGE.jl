function sympify_catch(exp::String,eqno::Int)
    try
        return Expr = Sym(exp)
    catch err
        if string(err.val) == "PyObject SympifyError()"
            # This is a terrible way to do this, but I can't figure a better way...
            println("Syntax in equation $(eqno) is invalid.")
            rethrow(err)
        else
            println("Error sympifying equation $(eqno):")
            rethrow(err)
        end
    end
end

function RegimeParameterDecomposition(parameter_values::Array{Float64,1},ergodic::Array{Float64,1})
    # outer constructor
    bar = dot(parameter_values,ergodic)
    hatarray = parameter_values - bar
    return RegimeParameterDecomposition(bar,hatarray)
end

function checksymbols(equations::Equations,vars::Variables,shocks::Shocks,parameters::Parameters)
    # This function checks all the expressions for symbols that are unrecognised.
    # Catching this early is helpful for avoiding errors later.
    validsymbols = vcat(vars.leads_sym,shocks.syms,vars.contemps_sym,vars.states_sym,parameters.generic_sym[:])
    for (i,expression) in enumerate(equations.syms)
        symbollist = SymPy.free_symbols(expression)
        for symbol in symbollist
            if !(symbol in validsymbols)
                error("Unknown symbol $(symbol) in equation $(i).")
            end
        end
    end
    return
end

function formatparametervalues(numregimes,numparameters,parameter_values,parameter_names)
    formatted = Array(Float64,numregimes,numparameters)
    switching = Array(Bool,length(parameter_names))
    fill!(switching,false)
    for (p,value,name) in zip(1:numparameters,parameter_values,parameter_names)
        if !isa(value,Number)
            switching[p] = true
            if length(value)!=numregimes
                error("Regime-switching parameter $(name) has the wrong number of declared parameters values. It must have one for each regime.")
            else
                formatted[:,p]=collect(value)
            end
        else
            formatted[:,p]=value
        end
    end
    return formatted,switching
end

function getergodic(transmatrix::Array{Float64,2})
    D,V = eig(transmatrix')
    eigindex = indmin(abs(D-1.0))
    ergodic = vec(V[:,eigindex]/sum(V[:,eigindex]))
    if any(abs(imag(ergodic)).>1e-16)
        error("Imaginary ergodic distribution")
    else
        return real(ergodic)
    end
end

function createtransmatrix_sym(numregimes::Int)
    out = Array(SymPy.Sym,numregimes)
    for j = 1:numregimes
        out[j] = Sym("TRANSITIONPROBABILITY_RX_RP$(j)")
    end
    return out
end

function sympifyallexpressions(expressions::Array{String,1},leads_sym::Array{SymPy.Sym,1},vars::Array{String,1})
    expressions_sym = Array(SymPy.Sym,length(expressions))

    lags = createsymversions(vars,"_lag") # Create a list of all lagged variable - these are _potential_ states
    writtenlags = createsymversions(vars,"(-1)") # Preallocating this saves time, by not having to repeatedly call Sym
    writtenleads = createsymversions(vars,"(+1)") # As above

    for (eqno,exp) in enumerate(expressions)
        Expr = sympify_catch(exp,eqno)
        expressions_sym[eqno] = Expr
    end
    
    subsargs = createpairedtuplelist([writtenlags writtenleads],[lags leads_sym])
    if !isempty(subsargs)
        expressions_sym = SymPy.subs(expressions_sym,subsargs...)
    end
    
    isstate = Array(Bool,length(lags))
    fill!(isstate,false)
    
    for expression in expressions_sym
        isstate = isstate | expressioncontains(expression,lags)
    end
    
    # Now create the states array and indices array
    numstates = countnz(isstate)
    states_sym = Array(SymPy.Sym,numstates)
    states_indices = Array(Int,numstates)
    counter = 0
    for (i,state) in enumerate(lags)
        if isstate[i]
            counter += 1
            states_sym[counter] = state
            states_indices[counter] = i
        end
    end

    return expressions_sym,isstate,states_sym,states_indices
end

function createsymversions(vars::Array{String,1},suffix=""::String)
    out = Array(SymPy.Sym,length(vars))
    for (i,var) in enumerate(vars)
        out[i] = Sym(string(var,suffix))
    end
    return out
end    

function createsymparameters(names::Array{String,1},switches::Array{Bool,1})
    out = Array(SymPy.Sym,2,length(names))
    for (i,name,switch) in zip(1:length(names),names,switches)
        if switch
            out[1,i] = Sym("$(name)_R")
            out[2,i] = Sym("$(name)_RP")
        else
            out[1,i] = Sym("$(name)")
            out[2,i] = out[1,i]
        end
    end
    return out
end

function evaluateSS(ssfunction::Array{Function,1},varvalues::Array{Float64,1},parametervalues::Array{Float64,1})
    out = Array(Float64,length(ssfunction))
    for (i,f) in enumerate(ssfunction)
        out[i] = f(varvalues...,parametervalues...)
    end
    return out
end

function lambdifyss(ss::Array{SymPy.Sym,1},vars::Array{SymPy.Sym,1},parameters::Array{SymPy.Sym,1})
    # Takes in sympy expressions, and symbolic object versions of the variables and parameters
    # Converts to a function that takes in arrays of variable values and parameter values
    # and returns equation errors
    functionarray = Array(Function,length(ss))
    for (i,expression) in enumerate(ss)
        functionarray[i] = SymPy.lambdify(expression,[vars;parameters])
    end
    return functionarray
end

function createSSparameters(numparameters::Int,parameter_affectsSS::Array{Bool,1},barversions::Array{SymPy.Sym},
			     parameter_names::Array{String},parameters::Array{SymPy.Sym},switching::Array{Bool,1})
    parametersSS = Array(SymPy.Sym,numparameters)
    for p = 1:numparameters
        if parameter_affectsSS[p]
            parametersSS[p] = barversions[p]
        elseif switching[p]
            # Define a new, non-state dependent symbol version
            parametersSS[p] = Sym(parameter_names[p])
        else
            # Just use the existing symbol
            parametersSS[p] = parameters[1,p]
        end
    end
    return parametersSS
end

function createSSsystem(numparameters,expressions,leads,lags,contemps,statesindices,shocks,
                        parameters,barversions,parameter_affectsSS,switching,parameter_names)

    parametersSS = createSSparameters(numparameters,parameter_affectsSS,barversions,parameter_names,parameters,switching)
    subsargs = createpairedtuplelist([leads; lags; shocks; reshape(parameters[1,:],numparameters); reshape(parameters[2,switching],countnz(switching))],
    [contemps; contemps[statesindices]; zeros(Float64,length(shocks)); parametersSS; parametersSS[switching]])
    expressionsSS = SymPy.subs(expressions,subsargs...)
    lambdifiedSS = lambdifyss(expressionsSS,contemps,parametersSS)
    return lambdifiedSS
end

function decomposeparameters(parameters,switching,affectsSS,parameter_values,ergodic,numparameters,numregimes)
    
    # Preallocate output arrays
    parameter_decomposition = Array(RegimeParameterDecomposition,numparameters)
    parameter_bar_sym = Array(SymPy.Sym,numparameters)
    parameter_regime_sym = Array(SymPy.Sym,numregimes,numparameters)
    
    sympyzero = Sym(0.0)
    
    for (p,name) in enumerate(parameters)
        if affectsSS[p]
            parameter_decomposition[p] = RegimeParameterDecomposition(parameter_values[:,p],ergodic)
            parameter_bar_sym[p] = Sym("$(name)_BAR")
            for r = 1:numregimes
                parameter_regime_sym[r,p] = Sym("$(name)_HAT_R$(r)")
            end
        else
            parameter_decomposition[p] = RegimeParameterDecomposition(0.0,zeros(Float64,1))
            parameter_bar_sym[p] = sympyzero
            if switching[p]
                for r = 1:numregimes
                    parameter_regime_sym[r,p] = Sym("$(name)_R$(r)")
                end
            else
                parameter_regime_sym[:,p] = Sym("$(name)")
            end
        end
    end
    
    return parameter_decomposition, parameter_bar_sym, parameter_regime_sym
end

function findparametersaffectSS(numparameters,parameter_names,parameter_sym,parameter_switching,
    parameter_values,expressions,leads,lags,contemps,statesindices,shocks,ss_guess,ergodic)
    
    println("Checking which of the switching parameters affect the steady state:")
    affectSS = Array(Bool,numparameters)
    fill!(affectSS,false)
    
    # generate bar versions for _all_ parameters
    baselinevalues = squeeze(permutedims(parameter_values,[2 1])*reshape(ergodic,(length(ergodic),1)),2)
    
    # Create and solve the SS using those baseline values
    SSsystem = createSSsystem(numparameters,expressions,leads,lags,contemps,statesindices,shocks,
    parameter_sym,Array(SymPy.Sym,0),affectSS,parameter_switching,parameter_names)

    sswrapper = x->evaluateSS(SSsystem,x,baselinevalues)

    print("Solving for steady state...")
    initialeval = sswrapper(ss_guess)
    if any(isnan(initialeval)) || any(initialeval.== Inf) || any(initialeval.== -Inf)
        error("Steady state guess does not evaluate (e.g. due to divide-by-zero). Supply a different initial guess.")
    end
    nlout = nlsolve(not_in_place(sswrapper), ss_guess)
    if !converged(nlout)
        println("FAILED.")
        println("Residuals at initial guess:")
        println(sswrapper(ss_guess))
        error("Could not find steady state.")
    end
    SSbaseline = nlout.zero
    println("done.")

    for (p,switches) in enumerate(parameter_switching)
        if switches
            print(string(parameter_names[p],"..."))
            values = parameter_values[:,p]
            # Find a value in the regime switching parameter values that is not equal to the baseline
            newindex = find(values.!=baselinevalues[p])
            if length(newindex) == 0
                println("does not.")
                warn("You declared $(parameter_names[p]) to be switching by giving it multiple values, but the values are all identical. This is very inefficient.")
            else
                newvalue = values[minimum(newindex)]
                newvalues = copy(baselinevalues)
                newvalues[p] = newvalue
                errors = evaluateSS(SSsystem,SSbaseline,newvalues)
                if any(abs(errors).>nlout.ftol)
                    affectSS[p] = true
                    println("does.")
                else
                    println("does not.")
                end
            end
        else
            affectSS[p] = false
        end
    end
    
    return affectSS,SSbaseline
end

function appearsmorethanonce(var,system)
    appeared = false
    for expression in system
        if var in SymPy.free_symbols(expression)
            if appeared
                return true
            else
                appeared = true
            end
        end
    end
    return false
end

function findstaticexpressions(expressions_sym,leads_sym,contemps_sym,states_sym,var_isstate,eq_isexo,var_isexo)
    numvars = length(contemps_sym)
    var_isstatic = copy(!var_isstate)
    expression_isstatic = trues(length(expressions_sym))
    for (i,expression) in enumerate(expressions_sym)
        for (j,lead) in enumerate(leads_sym)
            if lead in SymPy.free_symbols(expression)
                var_isstatic[j] = false
                expression_isstatic[i] = false
            end
        end
        # Only check this if haven't already flipped it
        if expression_isstatic[i]
            for state in states_sym
                if state in SymPy.free_symbols(expression)
                    expression_isstatic[i] = false
                end
            end
        end
    end
    
    # Finally, if a static variable appears only _once_ in the whole system, the equation it appears in has to be a static
    tocheck = Array(Bool,numvars)
    fill!(tocheck,false)
    for (i,var) in enumerate(contemps_sym)
        if var_isstatic[i] && !appearsmorethanonce(var,expressions_sym)
            tocheck[i] = true
        end
    end
    for var in contemps_sym[tocheck]
        for (i,expression) in enumerate(expressions_sym)
            if var in SymPy.free_symbols(expression)
                expression_isstatic[i] = true
            end
        end
    end
        
    # Now negate out any exos, since they'll be solved separately
    expression_isstatic = expression_isstatic & !eq_isexo
    var_isstatic = var_isstatic & !var_isexo

    return expression_isstatic,var_isstatic
end

function expressioncontains(expression,symbols,shortcircuit::Bool)
    # This is a short-circuit version of expression contains that is equivalent but quicker than any(expressioncontains(expression,symbols))
    for symbol in symbols
        if symbol in SymPy.free_symbols(expression)
            return true
        end
    end
    return false
end

function expressioncontains(expression,symbols)
    count = 0
    isin = Array(Bool,length(symbols))
    fill!(isin,false)
    for (i,symbol) in enumerate(symbols)
        if symbol in SymPy.free_symbols(expression)
            isin[i] = true
        end
    end
    return isin
end

function checkifcontainsforward(expressions,leads)
    contains = Array(Bool,length(expressions))
    fill!(contains,false)
    for (i,expression) in enumerate(expressions)
        contains[i] = expressioncontains(expression,leads,true)
    end
    return contains
end

function findexoequations(expressions,leads,contemps)
    contains_forward = checkifcontainsforward(expressions,leads)
    restrictedindices = collect(1:length(expressions))
    restrictedsystem = expressions[!contains_forward]
    restrictedindices = restrictedindices[!contains_forward]
    contempindices = collect(1:length(contemps))
    
    endocontemps = copy(contemps)
    exocontemps = Array(Int,0)
    exoequations = Array(Int,0)
    
    somethingchanged = true
    while somethingchanged
        somethingchanged = false
        keepexpression = trues(length(restrictedsystem))
        for (i,(expression,expression_index)) in enumerate(zip(restrictedsystem,restrictedindices))
            indices = expressioncontains(expression,endocontemps)
            if countnz(indices) == 1
                # Now figure out which one
                newexo = contempindices[indices]
                # Remove it from the list of endos
                endocontemps = endocontemps[!indices]
                contempindices = contempindices[!indices]
                push!(exocontemps,newexo[1])
                push!(exoequations,expression_index)
                keepexpression[i] = false
                somethingchanged = true
            end
        end
        # Now replace restricted system with just the equation we want to keep
        restrictedsystem = restrictedsystem[keepexpression]
        restrictedindices = restrictedindices[keepexpression]
    end
    
    # Convert to boolean arrays
    expression_isexo = Array(Bool,length(expressions))
    fill!(expression_isexo,false)
    expression_isexo[exoequations] = true
    var_isexo = Array(Bool,length(contemps))
    fill!(var_isexo,false)
    var_isexo[exocontemps] = true
    if countnz(expression_isexo) != countnz(var_isexo)
        error("Something has done wrong. Found $(countnz(expression_isexo)) exogenous equations, but $(countnz(var_isexo)) exogenous vars.")
    end
    
    return expression_isexo,var_isexo
end

function creategenerichats(names,affectsSS)
    generics = Array(SymPy.Sym,2,length(affectsSS))
    sympyzero = Sym(0.0)
    for (i,name) in enumerate(names)
        if affectsSS[i]
            generics[1,i] = Sym("$(name)_HAT_R")
            generics[2,i] = Sym("$(name)_HAT_RP")
        else
            generics[:,i] = sympyzero
        end
    end
    return generics
end

function Parameters!(meta,equations,vars,shocks,transmatrix,names,values,ss_guess)
    dictionary = Dict(zip(names,1:meta.numparameters))
    values,isswitching = formatparametervalues(meta.numregimes,meta.numparameters,values,names)
    generic_sym = createsymparameters(names,isswitching) 
    
    affectsSS,ss_values = findparametersaffectSS(meta.numparameters,names,generic_sym,isswitching,
    values,equations.syms,vars.leads_sym,vars.states_sym,vars.contemps_sym,vars.states_indices,shocks.syms,ss_guess,transmatrix.ergodic)
    
    generic_hat = creategenerichats(names,affectsSS)
    
    ss_guess[:] = ss_values # Overwrite the ss_guess with the actual solved for steady state
    # Parameter decomposition and create the sym versions
    decomposition,bar_sym,specificregime_sym = decomposeparameters(names,isswitching,affectsSS,values,transmatrix.ergodic,meta.numparameters,meta.numregimes)
    # Create a values matrix that subs in bars
    valueswithbar = copy(values)
    for p in 1:meta.numparameters
        if affectsSS[p]
            valueswithbar[:,p] = decomposition[p].bar
        end
    end
    
    return Parameters(names,values,dictionary,isswitching,affectsSS,decomposition,valueswithbar,generic_sym,generic_hat,bar_sym,specificregime_sym)
end

function TransitionMatrix(trans)
    # Sanity check the regime transition matrix
    numregimes = size(trans,1)
    if size(trans,2)!=size(trans,1)
        error("Regime transition matrix is not square.")
    end
    if any((abs(sum(trans,2))-1.0).>1e-12)
        error("Rows of the regime transition matrix do not sum to 1.")
    end
    # Create a sympy version of the regime trans matrix
    syms = createtransmatrix_sym(numregimes)
        
    # Get the ergodic distribution of the markovian regime process
    ergodic = getergodic(trans)
    return TransitionMatrix(trans,syms,ergodic,Sym("TRANSITIONPROBABILITY"))
end

function Shocks(names)
    dictionary = Dict(zip(names,1:length(names)))
    syms = createsymversions(names)
    return Shocks(names,dictionary,syms)
end

function rearrangeequations(equations)
    out = similar(equations)
    for (i,eq) in enumerate(equations)
        numequals = length(matchall(r"=", eq))
        if numequals > 1
            error("You have more than 1 equals sign in equation $(i).")
        elseif numequals == 1
            index = search(eq,'=')
            out[i] = string(eq[1:(index-1)]," - (", eq[(index+1):end] , ")")
        else
            out[i] = eq
        end
    end
    return out
end

function parseequationsandvars(vars,equations)
    
    if length(equations)!=length(vars)
        error("You have $(length(equations)) for $(length(vars)) endogenous variables.")
    end
    var_dictionary = Dict(zip(vars,1:length(vars)))
    # Define arrays for all the symbolic objects, except states, because we don't know them yet
    leads_sym = createsymversions(vars,"_lead")
    contemp_sym = createsymversions(vars)
    
    # Check the equations have zero or one equals signs. If one, convert to a one-sided expressions
    zeroedequations = rearrangeequations(equations)
    
    # Sympify the expressions and detect state variables while doing so
    eqsyms,var_isstate,states_sym,states_indices = sympifyallexpressions(zeroedequations,leads_sym,vars)
    
    # Split the equations into 'exogenous' and everything else. The exogenous ones are equations that depend entirely on
    # states and shocks - e.g. stochastic equations.
    eq_isexo,var_isexo = findexoequations(eqsyms,leads_sym,contemp_sym)
    # Now go through and mark which equations and variables are static- this gets used later to sub out non-transition equations
    eq_isstatic,var_isstatic = findstaticexpressions(eqsyms,leads_sym,contemp_sym,states_sym,var_isstate,eq_isexo,var_isexo)
    
    equations = Equations(equations,eq_isstatic,eq_isexo,eqsyms)
    vars = Variables(vars,var_dictionary,var_isstate,var_isstatic,var_isexo,states_indices,states_sym,leads_sym,contemp_sym)
    return equations,vars
end

function SteadyState(meta,equations,vars,parameters,shocks,ss_guess)
    if length(ss_guess) != meta.numvars
        error("Your supplied guess for the steady state has $(length(ss_guess)) entries, but you have $(meta.numvars) variables.")
    end
    ss_system = createSSsystem(meta.numparameters,equations.syms,vars.leads_sym,vars.states_sym,vars.contemps_sym,
    vars.states_indices,shocks.syms,parameters.generic_sym,parameters.bar_sym,parameters.affectsSS,parameters.isswitching,parameters.names)
    
    return SteadyState(ss_system,ss_guess)
end

function RSDSGEModel(var_names,shock_names,parameter_names,equations,parameter_values,regime_transmatrix,ss_guess)
    # Some useful constants first
    numvars = length(var_names)
    numshocks = length(shock_names)
    numparameters = length(parameter_names)
    numregimes = size(regime_transmatrix,1)
    
    # Check all the names are ok
    checknames(var_names)
    checknames(shock_names)
    checknames(parameter_names)
    
    transmatrix = TransitionMatrix(regime_transmatrix)
    shocks = Shocks(shock_names)
    equations,vars = parseequationsandvars(var_names,equations)
    meta = Meta(countnz(vars.isstate),numvars,numshocks,numregimes,numparameters,countnz(vars.isstatic),countnz(vars.isexo))
    parameters = Parameters!(meta,equations,vars,shocks,transmatrix,parameter_names,parameter_values,ss_guess)
    # Check for unknown symbols in any of the equations
    checksymbols(equations,vars,shocks,parameters)   

    # Create the steady state system
    steadystate = SteadyState(meta,equations,vars,parameters,shocks,ss_guess)

    # Create the perturbation system
    perturbationsystem = PerturbationSystem(meta,parameters,vars,shocks,equations,steadystate,transmatrix)
    
    return RSDSGEModel(meta,parameters,vars,shocks,equations,steadystate,transmatrix,perturbationsystem)
end

RSDSGEModel(v,s,p,e,pv,r) = RSDSGEModel(v,s,p,e,pv,r,zeros(length(v)))
