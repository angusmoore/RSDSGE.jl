function Base.copy(model::RSDSGEModel)
    return RSDSGEModel([deepcopy(getfield(model, k)) for k = 1:length(fieldnames(model))]...)
end

function Base.copy(meta::ModelMetaInfo)
    return Meta([deepcopy(getfield(meta, k)) for k = 1:length(fieldnames(meta))]...)
end

function Base.copy(parameters::Parameters)
    return Parameters([deepcopy(getfield(parameters, k)) for k = 1:length(fieldnames(parameters))]...)
end

function Base.copy(vars::Variables)
    return Variables([deepcopy(getfield(vars, k)) for k = 1:length(fieldnames(vars))]...)
end

function Base.copy(eqs::Equations)
    return Equations([deepcopy(getfield(eqs, k)) for k = 1:length(fieldnames(eqs))]...)
end

function Base.copy(ss::SteadyState)
    return SteadyState([deepcopy(getfield(ss, k)) for k = 1:length(fieldnames(ss))]...)
end

function Base.copy(trans::TransitionMatrix)
    return TransitionMatrix([deepcopy(getfield(trans, k)) for k = 1:length(fieldnames(trans))]...)
end

function Base.copy(pert::PerturbationSystem)
    return PerturbationSystem([deepcopy(getfield(pert, k)) for k = 1:length(fieldnames(pert))]...)
end

function checknames(names::Array{String})
    # This function throws an error for names that are unacceptable/reserved (i.e. pi, inf etc)
    badnames = ["Inf", "inf", "exp", "pi", "log", "beta"]

    for name in names
        if any(badnames.==name)
            error("Name: $(name) is a reserved name. Please change it.")
        end
    end
    return
end

function sympify_catch(exp::String,eqno::Int)
    try
        return Symbolics.parse_expr_to_symbolic(Meta.parse(replace(exp, " " => "")), EqEvalModule) # For some reason, spaces between brackets cause it to parse as tuple, rather than braces
    catch err
        @error("Syntax in equation $(eqno) is invalid.")
        rethrow(err)
    end
end

function checksymbols(equations::Equations,vars::Variables,shocks::Shocks,parameters::Parameters)
    # This function checks all the expressions for symbols that are unrecognised.
    # Catching this early is helpful for avoiding errors later.
    validsymbols = vcat(vars.leads_sym,shocks.syms,vars.contemps_sym,vars.states_sym,parameters.generic_sym[:])
    for (i,expression) in enumerate(equations.syms)
        symbollist = Symbolics.get_variables(expression)
        for symbol in symbollist
            if !(Symbol(symbol) in Symbol.(validsymbols))
                error("Unknown symbol $(symbol) in equation $(i).")
            end
        end
    end
    return
end
