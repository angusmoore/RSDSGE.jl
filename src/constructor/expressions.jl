
function sympifyallexpressions(expressions::Array{String,1},leads_sym::Array{SymPy.Sym,1},vars::Array{String,1})
    expressions_sym = Array(SymPy.Sym,length(expressions))

    lags = createsymversions(vars,"_lag") # Create a list of all lagged variable - these are _potential_ states
    writtenlags = createsymversions(vars,"(-1)") # Preallocating this saves time, by not having to repeatedly call Sym
    writtenleads = createsymversions(vars,"(+1)") # As above

    for (eqno,exp) in enumerate(expressions)
        Expr = sympify_catch(exp,eqno)
        expressions_sym[eqno] = Expr
    end
    
    subsargs = createpairedtuplelist([writtenlags writtenleads],[lags leads_sym])
    if !isempty(subsargs)
        expressions_sym = SymPy.subs(expressions_sym,subsargs...)
    end
    
    isstate = Array(Bool,length(lags))
    fill!(isstate,false)
    
    for expression in expressions_sym
        isstate = isstate | expressioncontains(expression,lags)
    end
    
    # Now create the states array and indices array
    numstates = countnz(isstate)
    states_sym = Array(SymPy.Sym,numstates)
    states_indices = Array(Int,numstates)
    counter = 0
    for (i,state) in enumerate(lags)
        if isstate[i]
            counter += 1
            states_sym[counter] = state
            states_indices[counter] = i
        end
    end

    return expressions_sym,isstate,states_sym,states_indices
end

function rearrangeequations(equations)
    out = similar(equations)
    for (i,eq) in enumerate(equations)
        numequals = length(matchall(r"=", eq))
        if numequals > 1
            error("You have more than 1 equals sign in equation $(i).")
        elseif numequals == 1
            index = search(eq,'=')
            out[i] = string(eq[1:(index-1)]," - (", eq[(index+1):end] , ")")
        else
            out[i] = eq
        end
    end
    return out
end
