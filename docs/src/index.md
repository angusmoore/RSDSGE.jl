# RSDSGE.jl Manual

```@meta
CurrentModule = RSDSGE
```

RSDSGE is a symbolic modelling engine for solving DSGE models. Itautomates the most tedious and error-prone part of solving DSGE models by perturbation methods: taking the perturbation. In that sense, it is very similar to existing tools such as Dynare. However, it is built with a more interactive workflow in mind than Dynare: models are ojbects, many of which can be kept in memory at once. This design allows you to experiment and test the parameters of your model interactively, at Julia's REPL, rather than executing a static 'model file'.

Models are solved one of two ways: the Grobner basis method of Foerster _et al_; or by numerical solver.

For a quick start guide, look at one of the exampls. The guide contains more detail on the commands used in the examples. The documentation lists and explains all of the functions that RSDSGE.jl makes available.

## Examples
```@contents
Pages = ["examples/FRWZ.md","examples/threeEQNK.md"]
Depth = 1
```

## Guide

```@contents
Pages = ["guide/models.md", "guide/dynare.md"]
```

## Documentation Index

```@index
Pages = ["lib/constructors.md", "lib/solve.md"]
Depth = 1
