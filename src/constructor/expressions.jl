
function sympifyallexpressions(expressions::Array{String,1},leads_sym::Array{Symbolics.Num,1},vars::Array{String,1})
    lags = createsymversions(vars,"_lag") # Create a list of all lagged variable - these are _potential_ states

    # Sub out the (-1) and (+1) syntax for _lag and _lead
    for var in vars
        expressions = map(x -> replace(x, var * "(+1)" => var * "_lead", var * "(-1)" => var * "_lag"), expressions)
    end

    expressions_sym = Symbolics.Num[sympify_catch(exp, eqno) for (eqno,exp) in enumerate(expressions)]

    isstate = fill(false, length(lags))

    for expression in expressions_sym
        isstate = isstate .| expressioncontains(expression,lags)
    end

    # Now create the states array and indices array
    numstates = count(!iszero, isstate)
    states_sym = Vector{Symbolics.Num}(undef, numstates)
    states_indices = Vector{Int}(undef, numstates)
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

function rearrangeequations(equations::Vector{String})
    out = similar(equations)
    for (i,eq) in enumerate(equations)
        numequals = length(findall(r"=", eq))
        if numequals > 1
            error("You have more than 1 equals sign in equation $(i).")
        elseif numequals == 1
            index = findfirst("=", eq)[1]
            out[i] = string(eq[1:(index-1)]," - (", eq[(index+1):end] , ")")
        else
            out[i] = eq
        end
    end
    return out
end
