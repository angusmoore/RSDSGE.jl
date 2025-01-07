function RegimeParameterDecomposition(parameter_values::Array{Float64,1},ergodic::Array{Float64,1})
    # outer constructor
    bar = LinearAlgebra.dot(parameter_values,ergodic)
    hatarray = parameter_values .- bar
    return RegimeParameterDecomposition(bar,hatarray)
end

function formatparametervalues(numregimes,numparameters,parameter_values,parameter_names)
    formatted = Array{Float64}(undef, numregimes,numparameters)
    switching = fill(false, length(parameter_names))
    for (p,value,name) in zip(1:numparameters,parameter_values,parameter_names)
        if !isa(value,Number)
            switching[p] = true
            if length(value)!=numregimes
                error("Regime-switching parameter $(name) has the wrong number of declared parameters values. It must have one for each regime.")
            else
                formatted[:,p]=collect(value)
            end
        else
            formatted[:,p] .= value
        end
    end
    return formatted,switching
end

function createsymversions(vars::Array{String,1},suffix=""::String)
    Symbolics.Num[Symbolics.variable(string(var,suffix)) for var in vars]
end

function createsymparameters(names::Array{String,1},switches::Vector{Bool})
    out = Array{Symbolics.Num}(undef, 2,length(names))
    for (i,name,switch) in zip(1:length(names),names,switches)
        if switch
            out[1,i] = Symbolics.variable("$(name)_R")
            out[2,i] = Symbolics.variable("$(name)_RP")
        else
            out[1,i] = Symbolics.variable("$(name)")
            out[2,i] = out[1,i]
        end
    end
    return out
end

function createSSparameters(numparameters::Int,parameter_affectsSS::Vector{Bool},barversions::Array{Symbolics.Num},
			     parameter_names::Array{String},parameters::Array{Symbolics.Num},switching::Vector{Bool})
    parametersSS = Vector{Symbolics.Num}(undef, numparameters)
    for p in 1:numparameters
        if parameter_affectsSS[p]
            parametersSS[p] = barversions[p]
        elseif switching[p]
            # Define a new, non-state dependent symbol version
            parametersSS[p] = Symbolics.variable(parameter_names[p])
        else
            # Just use the existing symbol
            parametersSS[p] = parameters[1,p]
        end
    end
    return parametersSS
end

function decomposeparameters(parameters,switching,affectsSS,parameter_values,ergodic,numparameters,numregimes)
    # Preallocate output arrays
    parameter_decomposition = Array{RegimeParameterDecomposition}(undef, numparameters)
    parameter_bar_sym = Array{Union{Symbolics.Num,Float64}}(undef, numparameters)
    parameter_regime_sym = Array{Union{Symbolics.Num,Float64}}(undef, numregimes,numparameters)

    for (p,name) in enumerate(parameters)
        if affectsSS[p]
            parameter_decomposition[p] = RegimeParameterDecomposition(parameter_values[:,p],ergodic)
            parameter_bar_sym[p] = Symbolics.variable("$(name)_BAR")
            for r = 1:numregimes
                parameter_regime_sym[r,p] = Symbolics.variable("$(name)_HAT_R$(r)")
            end
        else
            parameter_decomposition[p] = RegimeParameterDecomposition(0.0,zeros(Float64,1))
            parameter_bar_sym[p] = 0.0
            if switching[p]
                for r = 1:numregimes
                    parameter_regime_sym[r,p] = Symbolics.variable("$(name)_R$(r)")
                end
            else
                parameter_regime_sym[:,p] .= Symbolics.variable("$(name)")
            end
        end
    end

    return parameter_decomposition, parameter_bar_sym, parameter_regime_sym
end

function creategenerichats(names,affectsSS)
    generics = Array{Union{Symbolics.Num,Float64}}(undef, 2,length(affectsSS))
    for (i,name) in enumerate(names)
        if affectsSS[i]
            generics[1,i] = Symbolics.variable("$(name)_HAT_R")
            generics[2,i] = Symbolics.variable("$(name)_HAT_RP")
        else
            generics[:,i] .= 0.0
        end
    end
    return generics
end

function Parameters!(
    meta::ModelMetaInfo,
    equations::Equations,
    vars::Variables,
    shocks::Shocks,
    transmatrix::TransitionMatrix,
    names::Vector{String},
    values::Vector{Any},
    ss_guess::Vector{Float64}
)
    dictionary = Dict(zip(names,1:meta.numparameters))
    values,isswitching = formatparametervalues(meta.numregimes,meta.numparameters,values,names)
    generic_sym = createsymparameters(names,isswitching)

    affectsSS,ss_values = findparametersaffectSS(names,generic_sym,isswitching,
    values,equations.syms,vars.leads_sym,vars.states_sym,vars.contemps_sym,vars.states_indices,shocks.syms,ss_guess,transmatrix.ergodic)

    generic_hat = creategenerichats(names,affectsSS)

    ss_guess[:] = ss_values # Overwrite the ss_guess with the actual solved for steady state
    # Parameter decomposition and create the sym versions
    decomposition,bar_sym,specificregime_sym = decomposeparameters(names,isswitching,affectsSS,values,transmatrix.ergodic,meta.numparameters,meta.numregimes)
    # Create a values matrix that subs in bars
    valueswithbar = copy(values)
    for p in 1:meta.numparameters
        if affectsSS[p]
            valueswithbar[:,p] .= decomposition[p].bar
        end
    end

    return Parameters(names,values,dictionary,isswitching,affectsSS,decomposition,valueswithbar,generic_sym,generic_hat,bar_sym,specificregime_sym)
end
