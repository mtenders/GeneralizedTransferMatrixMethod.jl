![header](./docs/src/assets/banner.svg)

---

[![DOCS](https://img.shields.io/badge/docs-GeneralizedTransferMatrixMethod.jl-blue?style=flat-square)](https://mtenders.github.io/GeneralizedTransferMatrixMethod.jl/)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.7974658-blue?style=flat-square)](https://doi.org/10.5281/zenodo.7974658)



The present package implements the generalized transfer matrix
formalism as presented in [[1-3](#References)]. The authors provide a Python
implementation ([pyGTM](https://github.com/pyMatJ/pyGTM)), as well as a Matlab
implementation ([First version (2019)](https://doi.org/10.5281/zenodo.601496),
[Updated version (2020)](https://zenodo.org/record/3648041)).


## Installation

The package can be installed with the Julia package manager.
From the Julia REPL, type `]` to enter the Pkg REPL mode and run:

```
pkg> add GeneralizedTransferMatrixMethod
```

Or, equivalently, via the `Pkg` API:

```julia
julia> import Pkg; Pkg.add("GeneralizedTransferMatrixMethod")
```

## Example

This example of an glass–air interface is taken from the tutorial presented in the documentation

```julia
using GeneralizedTransferMatrixMethod
using LinearAlgebra
using Unitful: °, nm, μm, mm, m

# Define Materials
@permittivity "Glass" λ -> Diagonal(ones(3)) * 1.5
Air = Layer()

Interface = LayeredStructure(
    superstrate = Glass(),
    substrate = Air
)

# Calculate properties
λ = 1.55μm
α = 10°
ζ = sin(α) * sqrt(ϵ_Glass(λ)[1,1])

Properties = calculate_structure_properties(ζ, λ, Interface)

Rₚₚ, Rₛₛ, Rₚₛ, Rₛₚ = reflection(Properties)
Tₚ, Tₛ = transmission(ζ, Properties)
```

## References
1. [Passler, N. C. & Paarmann, A. Generalized 4 × 4 matrix formalism for light
   propagation in anisotropic stratified media: study of surface phonon
   polaritons in polar dielectric heterostructures. J. Opt. Soc. Am. B 34, 2128
   (2017)](http://doi.org/10.1364/JOSAB.34.002128). 
2. [Passler, N. C. & Paarmann, A. Generalized 4 × 4 matrix formalism for light
   propagation in anisotropic stratified media: study of surface phonon
   polaritons in polar dielectric heterostructures: erratum. J. Opt. Soc. Am. B
   36, 3246 (2019).](http://doi.org/10.1364/JOSAB.36.003246) 
3. [Passler, N. C., Jeannin, M. & Paarmann, A. Layer-resolved absorption of
   light in arbitrarily anisotropic heterostructures. Phys. Rev. B 101, 165425
   (2020).](http://doi.org/10.1103/PhysRevB.101.165425) 
