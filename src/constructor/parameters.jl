function RegimeParameterDecomposition(parameter_values::Array{Float64,1},ergodic::Array{Float64,1})
    # outer constructor
    bar = dot(parameter_values,ergodic)
    hatarray = parameter_values - bar
    return RegimeParameterDecomposition(bar,hatarray)
end

function formatparametervalues(numregimes::Int,numparameters::Int,parameter_values::Array,parameter_names::Array{String})
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

function decomposeparameters(parameters::Array{String,1},switching::Array{Bool,1},affectsSS::Array{Bool,1},parameter_values::Array{Float64,2},ergodic::Array{Float64,1},numparameters::Int,numregimes::Int)
    
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

function creategenerichats(names::Array{String,1},affectsSS::Array{Bool,1})
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

function Parameters!(meta::Meta,equations::Equations,vars::Variables,shocks::Shocks,transmatrix::TransitionMatrix,names::Array{String},values::Array,ss_guess::Array{Float64,1})
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
