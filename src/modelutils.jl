function findSS!(model::RSDSGEModel,guess::Array{Float64,1})
    ssparameters = [affects ? model.parameters.decomposition[p].bar : model.parameters.values[1, p] for (p,affects) in enumerate(model.parameters.affectsSS)]

    sswrapper = x->evaluateSS(model.steadystate.system,x,ssparameters)
    print("Solving for steady state...")
    nlout = NLsolve.nlsolve(NLsolve.not_in_place(sswrapper), guess)
    if !NLsolve.converged(nlout)
        println("FAILED.")
        println("Residuals at initial guess:")
        println(sswrapper(guess))
        error("Could not find steady state.")
    end
    model.steadystate.values = nlout.zero
    println("done.")
    return nlout.zero
end

function findSS!(model::RSDSGEModel)
    return findSS!(model,model.steadystate.values)
end

function updateparameters!(model::RSDSGEModel,name::String,value,reevaluateSS=true)
    index = model.parameters.dictionary[name]
    updateparameters!(model,index,value,reevaluateSS)
    return
end

function updateparameters!(model::RSDSGEModel,index::Int,value::Number,reevaluateSS=true)
    if model.parameters.isswitching[index]
        error("$(model.parameters.names[index]) is a state-dependent parameter, but you passed in only one value.")
    end
    model.parameters.values[:,index] .= value
    # Re-evaluate SS. I currently don't track whether non-switching parameters affect the SS or not. I should, because I could be much more
    # efficient if I don't need to re-evaluate the SS on all parameters
    if reevaluateSS
        findSS!(model)
    end
    return
end

function updateparameters!(model::RSDSGEModel,index::Int,value::Tuple,reevaluateSS=true)
    updateparameters!(model,index,collect(value),reevaluateSS)
    return
end

function updateparameters!(model::RSDSGEModel,index::Int,value::Array{Float64,1},reevaluateSS=true)
    if !model.parameters.isswitching[index]
        error("$(model.parameters.names[index]) is not a regime switching parameter, but you passed in many values for it.")
    elseif length(value) != model.meta.numregimes
        error("Dimension mismatch. You passed in $(length(value)) values for $(model.parameters.name[index]), but there are $(model.meta.numregimes) regimes.")
    end
    model.parameters.values[:,index] = value
    if model.parameters.affectsSS[index]
	# I need to re do the parameter decomposition here
	newdecomp = RegimeParameterDecomposition(model.parameters.values[:,index],model.transmatrix.ergodic)
	model.parameters.decomposition[index] = newdecomp
	if reevaluateSS
	    findSS!(model)
	end
    end
    return
end

function updateparameters!(model::RSDSGEModel,names::Array{String,1},values::Array{Any},reevaluateSS=true)
    if length(names) != length(values)
	error("You supplied a different number of parameters to change ($(length(names)) than values ($(length(values)).")
    end
    for (i,name,value) in zip(1:length(names),names,values)
	if i < length(names)
	    updateparameters!(model,name,value,false)
	else
	    updateparameters!(model,name,value,reevaluateSS)
	end
    end
    return
end
