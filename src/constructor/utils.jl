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

function createpairedtuplelist(keyarray::Array,valuearray::Array)
    out = ()
    for (key,val) in zip(keyarray,valuearray)
        out = (out...,(key,val))
    end
    return out
end

function checknames(names::Array{String})
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
