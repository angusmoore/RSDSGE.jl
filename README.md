# RSDSGE.jl

![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)
[![CI](https://github.com/angusmoore/RSDSGE.jl/workflows/CI/badge.svg)](https://github.com/angusmoore/RSDSGE.jl/actions/workflows/CI.yml/)
[![codecov](https://codecov.io/gh/angusmoore/RSDSGE.jl/graph/badge.svg?token=q19V8EAXY4)](https://codecov.io/gh/angusmoore/RSDSGE.jl)

A Julia package to solve regime switching dynamic stochastic general equilibrium (DSGE) models.

# Example

The model file `examples/RBC.txt` shows how to define the structure and parameters etc of the model:

```
vars: lammbda, r, Y, c, Inv, k, A
shocks: e_A

parameters: betta, sigma, delta, alpha, rho

parametervalues:
betta = 0.99
sigma = 1
delta = 0.025
alpha = 0.36
rho = (0.9,0.4)
end

equations:
exp(lammbda) = betta*exp(lammbda(+1))*(1+exp(r)-delta)
exp(r) = alpha*A*exp(k)^(alpha-1)
exp(Y) = exp(c) + exp(Inv)
exp(lammbda) = exp(c)^(-sigma)
exp(k) = exp(Inv) + (1-delta)*exp(k(-1))
log(A) = rho_R*log(A(-1)) + e_A
exp(Y) = A*exp(k)^alpha
end

transmatrix: [0.99 0.01
0.5 0.5]

ssguess:
lammbda = 1
c = 1
r = 1.01
k = 10
Inv = 1
Y = 2
A = 1
end
```

To load the model and solve it run:
```
using RSDSGE
model = RSDSGEModel(joinpath(pkgdir(RSDSGE), "examples", "RBC.txt"))
solve(model)
```