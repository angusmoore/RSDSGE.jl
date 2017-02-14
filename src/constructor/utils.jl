function copy(model::RSDSGEModel)
    return RSDSGEModel([deepcopy(getfield(model, k)) for k = 1:length(fieldnames(model))]...)
end

function copy(meta::Meta)
    return Meta([deepcopy(getfield(meta, k)) for k = 1:length(fieldnames(meta))]...)
end

function copy(parameters::Parameters)
    return Meta([deepcopy(getfield(parameters, k)) for k = 1:length(fieldnames(parameters))]...)
end

function copy(vars::Variables)
    return Meta([deepcopy(getfield(vars, k)) for k = 1:length(fieldnames(vars))]...)
end

function copy(eqs::Equations)
    return Meta([deepcopy(getfield(eqs, k)) for k = 1:length(fieldnames(eqs))]...)
end

function copy(ss::SteadyState)
    return Meta([deepcopy(getfield(ss, k)) for k = 1:length(fieldnames(ss))]...)
end

function copy(trans::TransitionMatrix)
    return Meta([deepcopy(getfield(trans, k)) for k = 1:length(fieldnames(trans))]...)
end

function copy(pert::PerturbationSystem)
    return Meta([deepcopy(getfield(pert, k)) for k = 1:length(fieldnames(pert))]...)
end

function createpairedtuplelist(keyarray::Array{Any},valuearray::Array{Any})
    out = ()
    for (key,val) in zip(keyarray,valuearray)
        out = (out...,(key,val))
    end
    return out
end

function checknames(names::Int)
    # This function throws an error for names that are unacceptable/reserved in sympy (i.e. pi, inf etc)
    badnames = Array(String,0)
    push!(badnames,"Inf")
    push!(badnames,"inf")
    push!(badnames,"exp")
    push!(badnames,"pi")
    push!(badnames,"log")
    push!(badnames,"beta")

    for name in names
        if any(badnames.==name)
            error("Name: $(name) is a reserved name in SymPy. Please change it.")
        end
    end
    return 
end
