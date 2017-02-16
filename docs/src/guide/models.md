# Dealing with Models

## The DSGEModel object

The core of the RSDSGE package is the RSDSGEModel type. This is a container for all of the information about your model: it's parameters, variables, shocks, equations and the system of perturbation equations.

Once constructed, a model is a variable like any other. You can copy it; have many in memory, corresponding to different things; update parameters on the fly; etc. This permits you to experiment with your model, interactively at the Julia REPL.

## Constructing a model

There are two ways to construct a model. The first uses a model file, which will be familiar to Dynare users. The file specifies parameters, parameter values, variables, shocks, the model's first order conditions and, optionally, a guess at the deterministic steady state. RSDGE takes it from there. This is a natural and easy way to write down your model and have the computer understand it.

The second way is to use the more direct constructor. Here you instead pass RSDSGE a series of arrays corresponding to parameters, parameter values etc. While somewhat less elegant than passing in a model file, it can be useful (for instance, constructing the model programatically where you have many sectors with the same first order conditions but different variable names).

The examples/ folder contains example model files. Model files are plain text files, and can have extension. To construct a model using a model file simply call (running this command requires your present working directory to be the package folder, try cd(Pkg.dir("RSDSGE"))):
```julia
RBC = DSGEModel("examples/RBC.txt")
'''


## Working with models



## Solving models

Solving a model is simple. Just call `solve(model)'. There are two methods available, using the keyword argument "method". Available options are "Grobner" (or, equivalently "FRWZ") to use the Grobner basis method from Foerster _et al_ or "numeric" to use a numeric solver. "Grobner" is the default and is recommended in most cases. As documented by Foerster _et al_, this method will find (and return) _all_ possible solutions to your model. A numeric solver will return only one, which may not be a sensible one. The Grobner basis method is typically more robust, although it can fail in large models.

`solve' returns a vector of RSDSGESolution. Each RSDSGE solution is a collection of matrices such that:
```julia
X_t = solution.MX*(X_{t-1} - SS) + solution.ME*E + solution.MC
'''
`SS' is the deterministic steady state, calculated as in Foerster _et al_ by taking the mean of regime-switching parameters under the ergodic distribution of the switching process. `E' is the column vector of shocks in the model.

Imaginary solutions are not returned. Each of the solutions is a candidate, and RSDSGE does not attempt to pick a `best' one for you. Foerster _et al_ discuss sensibe ways to choose among solutions on the basis of stability.

Higher-order solutions will be implemented in later versions.

As with model, you can have many solutions in memory. Once created, a solution corresponds to the model _as it was when you solved it_. If you subsequently update the model parameters, the solution does not change (you instead need to resolve the model).

This is useful for comparing solutions under different parameters, e.g.:
```julia
updateparameters!(model,"someparameter",0.9)
highparametersolution = solve(model)
updateparameters!(model,"someparameter",0.1)
lowparametersolution = solve(model)
'''

## IRFs and simulations

Yet to be implemented.

