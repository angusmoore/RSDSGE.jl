"""
   RSDSGEModel(file)
Construct a RSDSGE model from a model file. Pass in the path as a string to the file location.
For information on how to write a model file, see the model file documentation, or the examples contained in examples/
Returns your model as an instance of type RSDSGEModel.
"""
function RSDSGEModel(file::String)

    f = open(file)
	modelfilestring = readlines(f)
    close(f)

    # Find the indexes of important things
    ends = Bool[strip(line) == "end" for line in modelfilestring]
    ends = (1:length(modelfilestring))[ends]
    parameters = findline("parameters:",modelfilestring)
    parametervalues = findline("parametervalues:",modelfilestring)
    shocks = findline("shocks:",modelfilestring)
    vars = findline("vars:",modelfilestring)
    equations = findline("equations:",modelfilestring)
    transmatrix = findline("transmatrix: [",modelfilestring)
    transmatrixend = minimum((transmatrix:length(modelfilestring))[Bool[endswith(strip(line),"]") for line in modelfilestring[transmatrix:end]]])
    ssguess = findline("ssguess:",modelfilestring)

    equationend = minimum(ends[ends.>equations])
    ssguessend = minimum(ends[ends.>ssguess])
    valuesend = minimum(ends[ends.>parametervalues])

    vars = parsevars(modelfilestring[vars])
    shocks = parseshocks(modelfilestring[shocks])
    parameters = parseparameters(modelfilestring[parameters])
    parametervalues = parseparametervalues(modelfilestring[(parametervalues+1):(valuesend-1)],parameters)
    equations = parseequations(modelfilestring[(equations+1):(equationend-1)])
    transmatrix = parsetransmatrix(modelfilestring[transmatrix:transmatrixend])
    ssguess_dict = splitonequals(modelfilestring[(ssguess+1):(ssguessend-1)])

    ssguess = zeros(Float64,length(vars)) # Default to zero if unset
    for (key,val) in ssguess_dict
        idx = findfirst(x -> x == key, vars)
        if isnothing(idx)
            error("You have set a steady state guess for $key, but it is not a declared variable.")
        else
            ssguess[idx] = val
        end
    end

    return RSDSGEModel(string.(vars),string.(shocks),string.(parameters),string.(equations),parametervalues,transmatrix,ssguess)
end

function findline(keyword,modelfilestring)
    indices = Bool[startswith(strip(line),keyword) for line in modelfilestring]
    if count(!iszero, indices) != 1
	error("Model file contains $(count(!iszero, indices)) occurences of the $keyword keyword. It must appear only once.")
    end

    return (1:length(modelfilestring))[indices][1] # The extra one gets the scalar, rather than a 1 element array
end

function parsevars(vars)
    vars = strip(vars[6:end]) # Chop off the keyword and whitespace
    return splitoncomma(vars)
end

function parseshocks(shocks)
    shocks = strip(shocks[8:end]) # Chop off the keyword and whitespace
    return splitoncomma(shocks)
end

function parseparameters(parameters)
    parameters = strip(parameters[12:end]) # Chop off the keyword and whitespace
    return splitoncomma(parameters)
end

function parsetransmatrix(transmatrixstring)
    numregimes = length(transmatrixstring)
    transmatrix = Array{Float64}(undef,numregimes,numregimes)
    for r in 1:numregimes
	if r == 1
	    # First line, need to strip off the first bit
	    @warn "Make this more robust by finding [ instead of assuming it is 15"
	    transmatrixstring[r] = transmatrixstring[r][15:end]
	elseif r == numregimes
	    # Last line, strip off the closing ]
	    transmatrixstring[r] = strip(transmatrixstring[r])
	    transmatrixstring[r] = transmatrixstring[r][1:(end-1)]
	end
	vals = split(strip(transmatrixstring[r]), " ")
	for (rp,val) in enumerate(vals)
	    transmatrix[r,rp] = parse(Float64,val)
	end
    end
    return transmatrix
end

function parsetuple(input)
    if input[1] == '(' && input[end] == ')'
	input = input[2:(end-1)]
	vals = split(input,",")
	out = ()
	for val in vals
	    out = (out...,parse(Float64,val))
	end
	return out
    else
	println(input[1])
	println(input[end])
	error("$input is not valid tuple syntax.")
    end
end

function splitonequals(input)
    out = Dict{String,Union{Float64,Tuple}}()
    for i in eachindex(input)
	pair = split(input[i],"=")
	if length(pair)!=2
	    error("Error reading line '$(in[i])'; wrong number of equals signs.")
        else
	    if contains(pair[2],"(")
		out[strip(pair[1])] = parsetuple(strip(pair[2]))
	    else
		out[strip(pair[1])] = parse(Float64,strip(pair[2]))
	    end
        end
    end
    return out
end

function parseparametervalues(values,parameters)
    parametervalues = Array{Any}(undef, length(parameters))
    values_dict = splitonequals(values)
    unset = trues(length(parameters))
    for (key,val) in values_dict
        idx = findfirst(p -> p == key, parameters)
        if isnothing(idx)
            error("You have set a parameter value for $key, but it is not a declared parameter.")
        else
            parametervalues[idx] = val
	        unset[idx] = false
        end
    end
    if any(unset)
	    error("You have not set a value for parameter(s): $(parameters[unset]).")
    end
    return parametervalues
end

function parseequations(equations)
    map(strip, equations)
end


function splitoncomma(in)
    list = split(in,",")

    return map(strip, list)
end
