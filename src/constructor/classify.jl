function appearsmorethanonce(var::SymPy.Sym,system::Array{SymPy.Sym})
    appeared = false
    for expression in system
        if var in SymPy.free_symbols(expression)
            if appeared
                return true
            else
                appeared = true
            end
        end
    end
    return false
end

function findstaticexpressions(expressions_sym::Array{SymPy.Sym},leads_sym::Array{SymPy.Sym},contemps_sym::Array{SymPy.Sym},states_sym::Array{SymPy.Sym},var_isstate::BitArray,eq_isexo::BitArray,var_isexo::BitArray)
    numvars = length(contemps_sym)
    var_isstatic = copy(.!var_isstate)
    expression_isstatic = trues(length(expressions_sym))
    for (i,expression) in enumerate(expressions_sym)
        for (j,lead) in enumerate(leads_sym)
            if lead in SymPy.free_symbols(expression)
                var_isstatic[j] = false
                expression_isstatic[i] = false
            end
        end
        # Only check this if haven't already flipped it
        if expression_isstatic[i]
            for state in states_sym
                if state in SymPy.free_symbols(expression)
                    expression_isstatic[i] = false
                end
            end
        end
    end
    
    # Finally, if a static variable appears only _once_ in the whole system, the equation it appears in has to be a static
    tocheck = BitArray(numvars)
    fill!(tocheck,false)
    for (i,var) in enumerate(contemps_sym)
        if var_isstatic[i] && .!appearsmorethanonce(var,expressions_sym)
            tocheck[i] = true
        end
    end
    for var in contemps_sym[tocheck]
        for (i,expression) in enumerate(expressions_sym)
            if var in SymPy.free_symbols(expression)
                expression_isstatic[i] = true
            end
        end
    end
        
    # Now negate out any exos, since they'll be solved separately
    expression_isstatic = expression_isstatic .& .!eq_isexo
    var_isstatic = var_isstatic .& .!var_isexo

    return expression_isstatic,var_isstatic
end

function expressioncontains(expression::SymPy.Sym,symbols::Array{SymPy.Sym},shortcircuit::Bool)
    # This is a short-circuit version of expression contains that is equivalent but quicker than any(expressioncontains(expression,symbols))
    for symbol in symbols
        if symbol in SymPy.free_symbols(expression)
            return true
        end
    end
    return false
end

function expressioncontains(expression::SymPy.Sym,symbols::Array{SymPy.Sym})
    count = 0
    isin = BitArray(length(symbols))
    fill!(isin,false)
    for (i,symbol) in enumerate(symbols)
        if symbol in SymPy.free_symbols(expression)
            isin[i] = true
        end
    end
    return isin
end

function checkifcontainsforward(expressions::Array{SymPy.Sym},leads::Array{SymPy.Sym})
    contains = BitArray(length(expressions))
    fill!(contains,false)
    for (i,expression) in enumerate(expressions)
        contains[i] = expressioncontains(expression,leads,true)
    end
    return contains
end

function findexoequations(expressions::Array{SymPy.Sym},leads::Array{SymPy.Sym},contemps::Array{SymPy.Sym})
    contains_forward = checkifcontainsforward(expressions,leads)
    restrictedindices = collect(1:length(expressions))
    restrictedsystem = expressions[.!contains_forward]
    restrictedindices = restrictedindices[.!contains_forward]
    contempindices = collect(1:length(contemps))
    
    endocontemps = copy(contemps)
    exocontemps = Array{Int}(0)
    exoequations = Array{Int}(0)
    
    somethingchanged = true
    while somethingchanged
        somethingchanged = false
        keepexpression = trues(length(restrictedsystem))
        for (i,(expression,expression_index)) in enumerate(zip(restrictedsystem,restrictedindices))
            indices = expressioncontains(expression,endocontemps)
            if countnz(indices) == 1
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
    expression_isexo = BitArray(length(expressions))
    fill!(expression_isexo,false)
    expression_isexo[exoequations] = true
    var_isexo = BitArray(length(contemps))
    fill!(var_isexo,false)
    var_isexo[exocontemps] = true
    if countnz(expression_isexo) != countnz(var_isexo)
        error("Something has done wrong. Found $(countnz(expression_isexo)) exogenous equations, but $(countnz(var_isexo)) exogenous vars.")
    end
    
    return expression_isexo,var_isexo
end
