include("parameters.jl")
include("transmatrix.jl")
include("expressions.jl")
include("steadystate.jl")
include("classify.jl")




function Shocks(names)
    dictionary = Dict(zip(names,1:length(names)))
    syms = createsymversions(names)
    return Shocks(names,dictionary,syms)
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
