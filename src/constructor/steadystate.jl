
function evaluateSS(ssfunction::Array{Function,1},varvalues::Array{Float64,1},parametervalues::Array{Float64,1})
    map(f -> f(varvalues...,parametervalues...), ssfunction)
end

function lambdifyss(ss::Array{Symbolics.Num,1},vars::Array{Symbolics.Num,1},parameters::Array{Symbolics.Num,1})
    # Takes in sympy expressions, and symbolic object versions of the variables and parameters
    # Converts to a function that takes in arrays of variable values and parameter values
    # and returns equation errors
    return Function[Symbolics.build_function(expression,[vars;parameters]..., expression=Val{false},) for expression in ss]
end

function createSSsystem(numparameters::Integer,expressions::Vector{Symbolics.Num},leads::Vector{Symbolics.Num},lags::Vector{Symbolics.Num},contemps::Vector{Symbolics.Num},statesindices::Vector{<:Integer},shocks::Vector{Symbolics.Num},
                        parameters::Matrix{Symbolics.Num},barversions::Vector{Symbolics.Num},parameter_affectsSS::Vector{Bool},switching::Vector{Bool},parameter_names::Vector{String})

    parametersSS = createSSparameters(numparameters,parameter_affectsSS,barversions,parameter_names,parameters,switching)
    subsargs = Dict(
        old => new
        for (old, new) in zip(
            [leads; lags; shocks; reshape(parameters[1,:],numparameters); reshape(parameters[2,switching],count(!iszero, switching))],
            [contemps; contemps[statesindices]; zeros(Float64,length(shocks)); parametersSS; parametersSS[switching]]
        )
    )

    expressionsSS = map(x -> Symbolics.substitute(x, subsargs), expressions)

    lambdifiedSS = lambdifyss(expressionsSS,contemps,parametersSS)
    return lambdifiedSS
end

function findparametersaffectSS(
    parameter_names::Vector{String},
    parameter_sym::Matrix{Symbolics.Num},
    parameter_switching::Vector{Bool},
    parameter_values::Matrix{Float64},
    expressions::Vector{Symbolics.Num},
    leads::Vector{Symbolics.Num},
    lags::Vector{Symbolics.Num},
    contemps::Vector{Symbolics.Num},
    statesindices::Vector{Int},
    shocks::Vector{Symbolics.Num},
    ss_guess::Vector{Float64},
    ergodic::Vector{Float64}
)
    numparameters = length(parameter_names)

    println("Checking which of the switching parameters affect the steady state:")
    affectSS = fill(false, numparameters)

    # generate bar versions for _all_ parameters
    baselinevalues = reshape(permutedims(parameter_values,[2 1])*reshape(ergodic,(length(ergodic),1)), numparameters)

    # Create and solve the SS using those baseline values
    SSsystem = createSSsystem(numparameters,expressions,leads,lags,contemps,statesindices,shocks,
    parameter_sym,Vector{Symbolics.Num}(),affectSS,parameter_switching,parameter_names)

    sswrapper = x->evaluateSS(SSsystem,x,baselinevalues)

    print("Solving for steady state...")
    initialeval = sswrapper(ss_guess)
    if any(isnan.(initialeval)) || any(initialeval.== Inf) || any(initialeval.== -Inf)
        error("Steady state guess does not evaluate (e.g. due to divide-by-zero). Supply a different initial guess.")
    end
    nlout = NLsolve.nlsolve(NLsolve.not_in_place(sswrapper), ss_guess)
    if !NLsolve.converged(nlout)
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
            newindex = findfirst(x -> x != baselinevalues[p], values)
            if isnothing(newindex)
                println("does not.")
                @warn("You declared $(parameter_names[p]) to be switching by giving it multiple values, but the values are all identical. This is very inefficient.")
            else
                newvalue = values[newindex]
                newvalues = copy(baselinevalues)
                newvalues[p] = newvalue
                errors = evaluateSS(SSsystem,SSbaseline,newvalues)
                if any(abs.(errors).>nlout.ftol)
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

function SteadyState(meta::ModelMetaInfo,equations::Equations,vars::Variables,parameters::Parameters,shocks::Shocks,ss_guess::Vector{Float64})
    if length(ss_guess) != meta.numvars
        error("Your supplied guess for the steady state has $(length(ss_guess)) entries, but you have $(meta.numvars) variables.")
    end
    ss_system = createSSsystem(meta.numparameters,equations.syms,vars.leads_sym,vars.states_sym,vars.contemps_sym,
    vars.states_indices,shocks.syms,parameters.generic_sym,parameters.bar_sym,parameters.affectsSS,parameters.isswitching,parameters.names)

    return SteadyState(ss_system,ss_guess)
end
