function appearsmorethanonce(var::Symbolics.Num,system::Array{Symbolics.Num})
    appeared = false
    for expression in system
        if occursin(Symbolics.value(var), Symbolics.value(expression))
            if appeared
                return true
            else
                appeared = true
            end
        end
    end
    return false
end

function findstaticexpressions(expressions_sym::Array{Symbolics.Num},leads_sym::Array{Symbolics.Num},contemps_sym::Array{Symbolics.Num},states_sym::Array{Symbolics.Num},var_isstate::BitArray,eq_isexo::Vector{Bool},var_isexo::Vector{Bool})
    numvars = length(contemps_sym)
    var_isstatic = copy(.!var_isstate)
    expression_isstatic = trues(length(expressions_sym))
    for (i,expression) in enumerate(expressions_sym)
        for (j,lead) in enumerate(leads_sym)
            if occursin(Symbolics.value(lead), Symbolics.value(expression))
                var_isstatic[j] = false
                expression_isstatic[i] = false
            end
        end
        # Only check this if haven't already flipped it
        if expression_isstatic[i]
            for state in states_sym
                if occursin(Symbolics.value(state), Symbolics.value(expression))
                    expression_isstatic[i] = false
                end
            end
        end
    end

    # Finally, if a static variable appears only _once_ in the whole system, the equation it appears in has to be a static
    tocheck = fill(false, numvars)
    for (i,var) in enumerate(contemps_sym)
        if var_isstatic[i] && .!appearsmorethanonce(var,expressions_sym)
            tocheck[i] = true
        end
    end
    for var in contemps_sym[tocheck]
        for (i,expression) in enumerate(expressions_sym)
            if occursin(Symbolics.value(var), Symbolics.value(expression))
                expression_isstatic[i] = true
            end
        end
    end

    # Now negate out any exos, since they'll be solved separately
    expression_isstatic = expression_isstatic .& .!eq_isexo
    var_isstatic = var_isstatic .& .!var_isexo

    return expression_isstatic,var_isstatic
end

function expressioncontains(expression::Symbolics.Num,symbols::Array{Symbolics.Num},shortcircuit::Bool)
    # This is a short-circuit version of expression contains that is equivalent but quicker than any(expressioncontains(expression,symbols))
    for symbol in symbols
        if occursin(Symbolics.value(symbol), Symbolics.value(expression)) # Need to unwrap the Num (that's what value does)
            return true
        end
    end
    return false
end

function expressioncontains(expression::Symbolics.Num,symbols::Array{Symbolics.Num})
    Bool[occursin(Symbolics.value(symbol), Symbolics.value(expression)) for symbol in symbols] # The value is necessary to unwrap the Num
end

function checkifcontainsforward(expressions::Array{Symbolics.Num},leads::Array{Symbolics.Num})
    Bool[expressioncontains(expression,leads,true) for expression in expressions]
end

function findexoequations(expressions::Vector{Symbolics.Num},leads::Vector{Symbolics.Num},contemps::Array{Symbolics.Num})
    contains_forward = checkifcontainsforward(expressions,leads)
    restrictedindices = collect(1:length(expressions))
    restrictedsystem = expressions[.!contains_forward]
    restrictedindices = restrictedindices[.!contains_forward]
    contempindices = collect(1:length(contemps))

    endocontemps = copy(contemps)
    exocontemps = Vector{Int}()
    exoequations = Vector{Int}()

    somethingchanged = true
    while somethingchanged
        somethingchanged = false
        keepexpression = trues(length(restrictedsystem))
        for (i,(expression,expression_index)) in enumerate(zip(restrictedsystem,restrictedindices))
            indices = expressioncontains(expression,endocontemps)
            if count(!iszero, indices) == 1
                # Now figure out which one
                newexo = contempindices[indices]
                # Remove it from the list of endos
                endocontemps = endocontemps[.!indices]
                contempindices = contempindices[.!indices]
                push!(exocontemps,newexo[1])
                push!(exoequations,expression_index)
                keepexpression[i] = false
                somethingchanged = true
            end
        end
        # Now replace restricted system with just the equation we want to keep
        restrictedsystem = restrictedsystem[keepexpression]
        restrictedindices = restrictedindices[keepexpression]
    end

    # Convert to boolean arrays
    expression_isexo = fill(false, length(expressions))
    expression_isexo[exoequations] .= true
    var_isexo = fill(false, length(contemps))
    var_isexo[exocontemps] .= true
    if count(!iszero, expression_isexo) != count(!iszero, var_isexo)
        error("Something has done wrong. Found $(count(!iszero, expression_isexo)) exogenous equations, but $(count(!iszero, var_isexo)) exogenous vars.")
    end

    return expression_isexo,var_isexo
end
