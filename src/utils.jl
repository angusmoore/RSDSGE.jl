function copy(model::RSDSGEModel)
    return RSDSGEModel([deepcopy(getfield(model, k)) for k = 1:length(fieldnames(model))]...)
end

function copy(meta::Meta)
    return Meta([deepcopy(getfield(meta, k)) for k = 1:length(fieldnames(meta))]...)
end

function copy(parameters::Parameters)
    return Parameters([deepcopy(getfield(parameters, k)) for k = 1:length(fieldnames(parameters))]...)
end

function copy(vars::Variables)
    return Variables([deepcopy(getfield(vars, k)) for k = 1:length(fieldnames(vars))]...)
end

function copy(eqs::Equations)
    return Equations([deepcopy(getfield(eqs, k)) for k = 1:length(fieldnames(eqs))]...)
end

function copy(ss::SteadyState)
    return SteadyState([deepcopy(getfield(ss, k)) for k = 1:length(fieldnames(ss))]...)
end

function copy(trans::TransitionMatrix)
    return TransitionMatrix([deepcopy(getfield(trans, k)) for k = 1:length(fieldnames(trans))]...)
end

function copy(pert::PerturbationSystem)
    return PerturbationSystem([deepcopy(getfield(pert, k)) for k = 1:length(fieldnames(pert))]...)
end

function createpairedtuplelist(keyarray,valuearray)
    out = ()
    for (key,val) in zip(keyarray,valuearray)
        out = (out...,(key,val))
    end
    return out
end
