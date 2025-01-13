![header](./docs/src/assets/banner.svg)

---

[![DOCS](https://img.shields.io/badge/docs-GeneralizedTransferMatrixMethod.jl-blue?style=flat-square)](https://mtenders.github.io/GeneralizedTransferMatrixMethod.jl/)
[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.7974657-blue?style=flat-square)](https://doi.org/10.5281/zenodo.7974657)

This package implements the transfer matrix method  as presented in
[[1](#References)] for arbitrary relative permittivities ``\epsilon``,
*including non-reciprocal cases*. Arbitrary relative permeabilities ``\mu`` and
optical tensors ``\xi`` and ``\chi`` are also implemented but not yet tested.

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

Rₚₚ, Rₛₛ, Rₚₛ, Rₛₚ = calculate_reflection(λ, α, Interface)
Tₚₚ, Tₛₛ, Tₚₛ, Tₛₚ = calculate_transmission(λ, α, Interface)
```

## References
1. [Mackay, T. G. & Lakhtakia, A. The Transfer-Matrix Method in Electromagnetics and Optics. vol. 1 (2020).](https://doi.org/10.1007/978-3-031-02022-3)
