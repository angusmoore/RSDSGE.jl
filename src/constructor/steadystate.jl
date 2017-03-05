
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

function createSSsystem(numparameters,expressions,leads,lags,contemps,statesindices,shocks,
                        parameters,barversions,parameter_affectsSS,switching,parameter_names)

    parametersSS = createSSparameters(numparameters,parameter_affectsSS,barversions,parameter_names,parameters,switching)
    subsargs = createpairedtuplelist([leads; lags; shocks; reshape(parameters[1,:],numparameters); reshape(parameters[2,switching],countnz(switching))],
    [contemps; contemps[statesindices]; zeros(Float64,length(shocks)); parametersSS; parametersSS[switching]])
    expressionsSS = SymPy.subs(expressions,subsargs...)
    lambdifiedSS = lambdifyss(expressionsSS,contemps,parametersSS)
    return lambdifiedSS
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

function SteadyState(meta,equations,vars,parameters,shocks,ss_guess)
    if length(ss_guess) != meta.numvars
        error("Your supplied guess for the steady state has $(length(ss_guess)) entries, but you have $(meta.numvars) variables.")
    end
    ss_system = createSSsystem(meta.numparameters,equations.syms,vars.leads_sym,vars.states_sym,vars.contemps_sym,
    vars.states_indices,shocks.syms,parameters.generic_sym,parameters.bar_sym,parameters.affectsSS,parameters.isswitching,parameters.names)
    
    return SteadyState(ss_system,ss_guess)
end
